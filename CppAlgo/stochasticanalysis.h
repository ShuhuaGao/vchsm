#pragma once
#include <Eigen/dense>
#include "PicosElement.h"
#include <vector>
#include "types.h"
#include "StructDefinitions.h"


/**
 * Implementation of stochasticanalysis.m
 *@note picos is passed by reference which will be changed in this functions
*/
void stochasticanalysis(Eigen::Ref<const Eigen::TRowVectorX> x, Eigen::TFloat fs, int N, PicosStructArray& picos, int ordenLPC);
