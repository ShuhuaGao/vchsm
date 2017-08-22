#include <cstdio>
#include "train_C.h"
#include "convert_C.h"

int main()
{
	// train a model with the given 20 source and target speaker's audios
	// the audios files are named from 1 to 20
	const char* sourceAudioDir = "E:/GitHub/vchsm/Audios/source_train/";
	const char*  targetAudioDir = "E:/GitHub/vchsm/Audios/target_train/";
	const int numTrainSamples = 20;
	const char* sourceAudioList[numTrainSamples];
	const char* targetAudioList[numTrainSamples];
	for (int i = 0; i < numTrainSamples; ++i)
	{
		char* buff = new char[100];
		std::sprintf(buff, "%s%d.wav", sourceAudioDir, i + 1);
		sourceAudioList[i] = buff;
		buff = new char[100];
		std::sprintf(buff, "%s%d.wav", targetAudioDir, i + 1);
		targetAudioList[i] = buff;
	}
	// model file to be generated
	const char* modelFile = "E:/GitHub/vchsm/models/Model.dat";
	// start training	
	trainHSMModel(sourceAudioList, targetAudioList, numTrainSamples, 4, modelFile);
	// deallocate
	for (int i = 0; i < numTrainSamples; ++i)
	{
		delete[] sourceAudioList[i];
		delete[] targetAudioList[i];
	}
	// perform conversion
	const char* testAudio = "E:/GitHub/vchsm/Audios/test/jal_in_41_3.wav";
	const char* testAudioConverted = "E:/GitHub/vchsm/Audios/test/jal_in_41_3_c.wav";
	convertSingle(modelFile, testAudio, testAudioConverted);
	// now we can compare the above audio before and after conversion
}