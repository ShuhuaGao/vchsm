#pragma once
#include <vector>
#include "PicosElement.h"
#include "DynamicTimeWarping.h"
#include "StructDefinitions.h"

/**
 *Implementation of PrepareParallelData.
 *@param sourceHSMFeatureList the picos list for all the source training samples
 *@param targetHSMFeatureList the picos list for all the target training samples
 *@return p1v and p2v in PrepareParallelData.m
 *@note Since we don't use .mat files for intermediate data exchange, we just pass these data as function arguments.
*/
std::pair<PicosStructArray, PicosStructArray> PrepareParallelData(const std::vector<PicosStructArray>& sourceHSMFeatureList, const std::vector<PicosStructArray>& targetHSMFeatureList, int index);
