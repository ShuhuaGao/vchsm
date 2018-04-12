#pragma once

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
    /**
     * Train the HSM model with given audios and write the model parameters into given text files
     *@param configFile a configuration text file to specify the training process
     * The configuration file should presented in the following form
     *---------------start of the config file-----------------
     * number_of_training_samples number_of_GMM
     * source_audio_1 target_audio_1
     * source_audio_2 target_audio_2
     *[...more source and target audios here....]
     * model_file_path_to_store_the_trained_parameters
     *---------------end of the config file----------------------
	 *@param verbose: set a non-zero integer if you want more information displayed during training
     *@note The input source and target audios must be: (1) .wav type; (2) 16bits per sample, 16KHz sampling frequency and a single channel
     */
    void trainHSMModelByConfig(const char* configFile, int verbose);
    
    
    /**
     * Train an HSM model
     *@param sourceAudioList the training .wav files of the source speaker
     *@param targetAudioList the training .wav files of the target speaker
     *@param numTrainSamples number of the training samples (aka, the size of sourceWavFiles and targetWavFiles)
     *@param m number of Gaussian components
     *@param modelFile path of the model file to be created
	 *@param verbose: set a non-zero integer if you want more information displayed during training
     *@note The input source and target audios must be: (1) .wav type; (2) 16bits per sample; (3) 16KHz sampling frequency and (4) a single channel (mono).
     */
    void trainHSMModel(const char* sourceAudioList[], const char* targetAudioList[], 
		int numTrainSamples, int m, const char* modelFile, int verbose);

#ifdef __cplusplus
}
#endif
