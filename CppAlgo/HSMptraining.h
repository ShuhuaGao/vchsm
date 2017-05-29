#pragma once
#include "PicosElement.h"
#include "StructDefinitions.h"
#include <vector>
#include <utility>

/**
 * Implementation of HSMptraining.m
 *@param corpus corresponding to the corpus.mat in MATLAB, including of two variables: p1v and p2v. Please check PrepareParallelData.m.
 *@param m number of gaussian components of the function (def=8)
 *@return the HSM model parameters, including the forward and the inverse one.
 *@note We don't use .mat files for data exchange here, therefore the 2nd parameter in HSMptraining.m is neglected.
 *@remark In the conversion part, the inverse model is unused. No need to generate it.
*/
HSMModel HSMptraining(const std::pair<PicosStructArray, PicosStructArray>& corpus, int m);
