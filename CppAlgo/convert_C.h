#pragma once



#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

	/**
	* Convert a single wav file.
	*@param modelFile specify the binary file storing the training parameters
	*@param wavFile specify the original wav file to be converted
	*@param convertedWavFile specify path for the converted wav file
	*@note This algorithm only supports wav audio files with single channel, 16 bits and 16KHz sampling frequency.
	*/
	void convertSingle(const char* modelFile, const char* wavFile, const char* convertedWavFile);

	/**
	 * Convert multiple wav audios in batch mode using a configuration file.
	 *@param conversionConfigFile a text file to configure the conversion process
	 *@note The conversion configuration file should be presented in the following form:
	 *	number_of_source_audio_files_to_be_converted
	 *	source_audio_file_1 converted_audio_file_1
	 *	source_audio_file_2 converted_audio_file_2
	 *	source_audio_file_3 converted_audio_file_3
	 *	[more source and converted audio files]
	 *	model_file_which_contains_training_parameters
	*/
	void convertBatchByConfig(const char* conversionConfigFile);
    
    /**
     * Convert multiple wav audios in batch mode
     *@param modelFile specify the binary file storing the training parameters
     *@param audioFileList specify the original wav files to be converted
     *@param convertedAudioFileList specify path for the converted wav files
     *@param numAudios the total number of audios for conversion, i.e., the size of wavFileList and convertedWavFileList
     *@note This algorithm only supports wav audio files with single channel, 16 bits and 16KHz sampling frequency.
     */
    void convertBatch(const char* modelFile, const char* audioFileList[], const char* convertedAudioFileList[], int numAudios);


#ifdef __cplusplus
}
#endif // __cplusplus
