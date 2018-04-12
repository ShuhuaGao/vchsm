#pragma once
#include <Eigen/Dense>
#include <tuple>
#include "types.h"

/**
 * Implementation of mfcc function in mfcc.m
 *@note Here there is no const qualifier for speech intentionally since it may be changed
 *In mfcc.m, there are three returns. However, only the first is ever used.
*/
Eigen::TMatrixX mfcc(Eigen::Ref<Eigen::TVectorX> speech, Eigen::TFloat fs);
