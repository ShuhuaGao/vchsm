#include "train_C.h"
#include "types.h"
#include "WavRW.h"
#include "StructDefinitions.h"
#include "HSManalyze_mfcc.h"
#include "HSMptraining.h"
#include "PrepareParallelData.h"
#include "modelSerialization.h"
#include <vector>
#include <fstream>
#include <algorithm>
#include <dispatch/dispatch.h>
#ifdef SHOWINFO
#include <iostream>
#endif

// Use Apple GCD parallel framework to accelarate the training process, see http://aras-p.info/blog/2013/02/03/parallel-for-in-apples-gcd/

// the number of hardware cores in the CPU
#define CPU_CORE_NUMBER 4

void trainHSMModelByConfig(const char* configFile)
{

	// 0. parse the training configuration file
	std::ifstream inFile(configFile);
	int numTrainSamples, m;
	inFile >> numTrainSamples >> m;
	std::vector<std::string> sourceAudioFiles(numTrainSamples), targetAudioFiles(numTrainSamples);
	for (int i = 0; i < numTrainSamples; i++)
	{
		inFile >> sourceAudioFiles[i] >> targetAudioFiles[i];
	}
	std::string modelFile;
    inFile >> modelFile;
	inFile.close();
    
    // 1. train
    std::vector<const char*> sourceAudioList(sourceAudioFiles.size());
    std::transform(sourceAudioFiles.cbegin(), sourceAudioFiles.cend(), sourceAudioList.begin(),
                   [](const std::string& file){return file.c_str();});
    std::vector<const char*> targetAudioList(targetAudioFiles.size());
    std::transform(targetAudioFiles.cbegin(), targetAudioFiles.cend(), targetAudioList.begin(),
                   [](const std::string& file){return file.c_str();});
    trainHSMModel(sourceAudioList.data(), targetAudioList.data(), numTrainSamples, m, modelFile.c_str());
}

void trainHSMModel(const char* sourceAudioList[], const char* targetAudioList[], int numTrainSamples, int m, const char* modelFile)
{
#ifdef SHOWINFO
    std::cout << "Using SIMD:\n" << Eigen::SimdInstructionSetsInUse() << std::endl << std::endl;
#endif

    // 1. read the audio
#ifdef SHOWINFO
    std::cout << "Read audios for training..." << std::endl;
#endif
    __block std::vector<Eigen::TRowVectorX> sourceAudioDataList, targetAudioDataList;
    sourceAudioDataList.reserve(numTrainSamples);
    targetAudioDataList.reserve(numTrainSamples);
    for (int i = 0; i < numTrainSamples; i++)
    {
        sourceAudioDataList.emplace_back(readWav(sourceAudioList[i]));
        targetAudioDataList.emplace_back(readWav(targetAudioList[i]));
    }
    
    // 2. MFCC
    __block std::vector<PicosStructArray> sourceHSMFeatureList(numTrainSamples), targetHSMFeatureList(numTrainSamples);
    const int fs = 16000;
    // dispatch in parallel MFCC
    int gcdWorkStride = numTrainSamples / int(CPU_CORE_NUMBER);
    dispatch_apply(CPU_CORE_NUMBER, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^(std::size_t strideIndex){
        int start = (int)strideIndex * gcdWorkStride;
        for (int i = start; i < start + gcdWorkStride; i++)
        {
#ifdef SHOWINFO
            std::cout << "MFCC: " << sourceAudioList[i] << std::endl;
#endif
            sourceHSMFeatureList[i] = HSManalyze_mfcc(sourceAudioDataList[i], fs);
#ifdef SHOWINFO
            std::cout << "MFCC: ";
            std::cout << targetAudioList[i] << std::endl;
#endif
            targetHSMFeatureList[i] = HSManalyze_mfcc(targetAudioDataList[i], fs);
        }
    });
    // handle the possibly remaining samples in serial mode
    int start = int(CPU_CORE_NUMBER) * gcdWorkStride;
    for (int i = start; i < numTrainSamples; i++)
    {
#ifdef SHOWINFO
        std::cout << "MFCC: " << sourceAudioList[i] << std::endl;
#endif
        sourceHSMFeatureList[i] = HSManalyze_mfcc(sourceAudioDataList[i], fs);
#ifdef SHOWINFO
        std::cout << "MFCC: ";
        std::cout << targetAudioList[i] << std::endl;
#endif
        targetHSMFeatureList[i] = HSManalyze_mfcc(targetAudioDataList[i], fs);
    }
    
    //3. prepare parallel data with the above HSM features to get the parallel corpus
#ifdef SHOWINFO
    std::cout << "Prepare parallel..." << std::endl;
#endif
    auto corpus = PrepareParallelData(sourceHSMFeatureList, targetHSMFeatureList, numTrainSamples);
    
    // 4. use the corpus data to perform training
#ifdef SHOWINFO
    std::cout << "Training..." << std::endl;
#endif
    auto model = HSMptraining(corpus, m);
    
    // 5. save model into text files
#ifdef SHOWINFO
    std::cout << "Writing model file..." << std::endl;
#endif
    serializeModel(model, modelFile);
}
