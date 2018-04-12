#include "../include/vchsm/train_C.h"
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
#include <iostream>

void trainHSMModelByConfig(const char* configFile, int verbose)
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
    trainHSMModel(sourceAudioList.data(), targetAudioList.data(), numTrainSamples, m, modelFile.c_str(), verbose);
	

}

void trainHSMModel(const char* sourceAudioList[], const char* targetAudioList[], 
	int numTrainSamples, int m, const char* modelFile, int verbose)
{
	if (verbose != 0)
	{
		std::cout << ">> Train...\n";
		std::cout << "Using SIMD: " << Eigen::SimdInstructionSetsInUse() << std::endl;
	}

    // 1. read the audio
	if (verbose != 0)
		std::cout << "\t[1] Read audios for training..." << std::endl;
    std::vector<Eigen::TRowVectorX> sourceAudioDataList, targetAudioDataList;
    sourceAudioDataList.reserve(numTrainSamples);
    targetAudioDataList.reserve(numTrainSamples);
    for (int i = 0; i < numTrainSamples; i++)
    {
		if (verbose != 0)
		{
			std::cout << "\t\t Read " << sourceAudioList[i] << std::endl;
			std::cout << "\t\t Read " << targetAudioList[i] << std::endl;
		}

        sourceAudioDataList.emplace_back(readWav(sourceAudioList[i]));
        targetAudioDataList.emplace_back(readWav(targetAudioList[i]));
    }
    
    // 2. MFCC
    std::vector<PicosStructArray> sourceHSMFeatureList(numTrainSamples), targetHSMFeatureList(numTrainSamples);
    const int fs = 16000;
	if (verbose != 0)
		std::cout << "\t[2] MFCC...\n";
    
    for (int i = 0; i < numTrainSamples; i++)
    {
		if (verbose != 0)
			std::cout << "\t\t MFCC: " << sourceAudioList[i] << std::endl;
        sourceHSMFeatureList[i] = HSManalyze_mfcc(sourceAudioDataList[i], fs);
		if (verbose != 0)
			std::cout << "\t\t MFCC: " << targetAudioList[i] << std::endl;
        targetHSMFeatureList[i] = HSManalyze_mfcc(targetAudioDataList[i], fs);
    }
    
    //3. prepare parallel data with the above HSM features to get the parallel corpus
	if (verbose != 0)
		std::cout << "\t[3] Prepare parallel data\n";
    auto corpus = PrepareParallelData(sourceHSMFeatureList, targetHSMFeatureList, numTrainSamples);
    
    // 4. use the corpus data to perform training
	if (verbose != 0)
		std::cout << "\t[4] Perform HSM training..." << std::endl;
    auto model = HSMptraining(corpus, m);
    
    // 5. save model into text files
	if (verbose != 0)
		std::cout << "\t[5] Writing model file..." << std::endl;
    serializeModel(model, modelFile);

	if (verbose != 0)
		std::cout << ">> All training finished. The produced model file is " << modelFile << ".\n";
}
