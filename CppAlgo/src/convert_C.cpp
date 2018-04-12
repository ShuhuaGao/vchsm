#include "../include/vchsm/convert_C.h"
#include "WavRW.h"
#include "modelSerialization.h"
#include "HSManalyze.h"
#include "HSMwfwconvert.h"
#include "HSMsynthesize.h"
#include <fstream>
#include <vector>
#include <iostream>


void convertSingle(const char * modelFile, const char * wavFile, const char * convertedWavFile, int verbose)
{
	if (verbose != 0) 
		std::cout << ">> Converting " << wavFile << "...\n";
	// 1. read the audio file into a vector
	auto x = readWav(wavFile);
	int L = static_cast<int>(x.size());

	// 2. read the trained model
	auto model = deserializeModel(modelFile);

	// 3. analyze
	int fs = 16000;
	auto picos = HSManalyze(x, fs);

	// 4. convert
	HSMwfwconvert(model, picos);

	// 5. synthesize
	auto y = HSMsynthesize(picos, L);

	// 6. write the converted data into a new audio file
	writeWav(y, convertedWavFile);
	if (verbose != 0)
		std::cout << ">> Finished. The converted audio file is " << convertedWavFile << ".\n";
}


void convertBatchByConfig(const char * conversionConfigFile, int verbose)
{
	// 1. parse the config file
	int numAudios = 0;
	std::vector<std::string> srcAudioList;
	std::vector<std::string> convertedAudioList;
	std::string modelFile;

	std::ifstream is;
	is.open(conversionConfigFile);
	is >> numAudios;
	std::string srcAudio, convertedAudio;
	for (int i = 0; i < numAudios; i++)
	{
		is >> srcAudio >> convertedAudio;
		srcAudioList.emplace_back(std::move(srcAudio));
		convertedAudioList.emplace_back(std::move(convertedAudio));
	}
	is >> modelFile;
	is.close();

	// 2. perform multiple conversions in parallel mode
	for (int i = 0; i < numAudios; i++)
	{
		convertSingle(modelFile.c_str(), srcAudioList[i].c_str(), convertedAudioList[i].c_str(), verbose);
	}

}

void convertBatch(const char* modelFile, const char* audioFileList[], const char* convertedAudioFileList[], 
	int numAudios, int verbose)
{
    for (int i = 0; i < numAudios; i++)
    {
        convertSingle(modelFile, audioFileList[i], convertedAudioFileList[i], verbose);
    }
}
