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
     *@note The input source and target audios must be: (1) .wav type; (2) 16bits per sample, 16KHz sampling frequency and a single channel
     */
    void trainHSMModelByConfig(const char* configFile);
    
    
    /**
     * Train an HSM model
     *@sourceAudioList the training .wav files of the source speaker
     *@targetAudioList the training .wav files of the target speaker
     *@numTrainSamples number of the training samples (aka, the size of sourceWavFiles and targetWavFiles)
     *@m number of Gaussian components
     *@modelFile path of the model file to be created
     *@note The input source and target audios must be: (1) .wav type; (2) 16bits per sample, 16KHz sampling frequency and a single channel
     */
    void trainHSMModel(const char* sourceAudioList[], const char* targetAudioList[], int numTrainSamples, int m, const char* modelFile);

#ifdef __cplusplus
}
#endif
