#pragma once
#include "StructDefinitions.h"

/**
 * Implmementation of HSMsynthesize.m
 *@param audioFeature feature of the converted audio to be synthesized 
 *@param L size of the original audio
 *@return a row vector of length L, which gives the data of the synthesized converted audio
*/
Eigen::TRowVectorX HSMsynthesize(const PicosStructArray& audioFeature, int L);