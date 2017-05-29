#pragma once
#include "types.h"

/**
 * Implementation of aaalsf.m.
 *@remark Since only the first return value in MATLAB is ever used, we only return the 1st one here.
*/
Eigen::TVectorX aaalsf(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::TFloat f0, int p = 14);
