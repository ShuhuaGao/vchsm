#pragma once
#include "types.h"
#include <string>
/**
 *@file Functions for reading and write .WAV files, whose sampling rate is 16KHz with 16bits per sample and a single channel.
 *@note Refer to http://stackoverflow.com/questions/13660777/c-reading-the-data-part-of-a-wav-file
*/

/**
* Read an audio file and gives its data in a vector form
*@param wavFile the audio file path
*@return the normalized audio data between [-1, 1]
*@note The audio must be 16bit, 16000KHz sampling and only has a single channel.
*/
Eigen::TRowVectorX readWav(const std::string & wavFile);

/**
 * Write the data into a WAV audio file
 *@param data the normalized audio data in [-1, 1]
 *@note The created audio file will be 16bit, 16000KHz sampling and only has a single channel.
*/
void writeWav(Eigen::Ref<const Eigen::TRowVectorX> data, const std::string& wavFile);

