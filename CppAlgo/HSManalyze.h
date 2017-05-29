#pragma once
#include "StructDefinitions.h"

/**
 * Implementation of HSManalyze.m
 *@param x data of the audio
 *@param fs the sampling frequency of the audio
 *@note There are some minor discriminations between HSManalyze.m and HSManalyze_mfcc.m.
*/
PicosStructArray HSManalyze(Eigen::Ref<Eigen::TRowVectorX> x, Eigen::TFloat fs);
