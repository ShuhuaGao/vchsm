#pragma once
#include "PicosElement.h"
#include <vector>
#include <Eigen/Dense>
#include "types.h"
#include "StructDefinitions.h"

/**
 * Implementation of HSManalyze_mfcc in HSManalyze_mfcc.m.
 *@param x data of the audio, normalized between -1 and 1
 *@param fs sampling frequency of the audio, which can only be 16000
 *@return an array of structures, each containing the model information
 *@note As specified in the HSManalyze_mfcc.m, only wav files sampled at 16KHz 16bits mono (single channel), are accepted.
 *@remark Unlike the original HSManalyze_mfcc function in MATLAB, we think it is better to separate the file I/O and data processing.
 *The role of this algorithm is to process input data and return the processed data, excluding the file I/O operations, since file I/O
 *is platform dependent. 
 *In HSManalyze_mfcc.m, the audio is first read and then the algorithm yields three variables: L, fs and picos, which is then saved into a mat file.
 *Here, this function returns the variable picos, which corresponds a mat file.
*/
PicosStructArray HSManalyze_mfcc(Eigen::Ref<Eigen::TRowVectorX> x, Eigen::TFloat fs);
