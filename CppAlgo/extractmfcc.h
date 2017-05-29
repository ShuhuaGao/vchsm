#pragma once
#include <Eigen/Dense>
#include <vector>
#include "PicosElement.h"
#include "types.h"
#include "StructDefinitions.h"

/**
 * Implementation of extractmfcc in extractmfcc.m.
*/
void extractmff(Eigen::Ref<Eigen::TRowVectorX> x, Eigen::TFloat fs, PicosStructArray& picos);
