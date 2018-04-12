#pragma once
#include <Eigen/Dense>
#include <tuple>
#include "types.h"

/**
 * Implementation of the trifbank function in trifbank.m.
 *@note The last two parameters (function handle type in MATLAB) are neglected here since this function is only called in mfcc.m.
 *And the last two parameters are just two fixed functions.
 *It is further found that only the first return of trifbank.m is used.
*/
Eigen::TMatrixX trifbank(int M, int K, const Eigen::TRowVector2& R, Eigen::TFloat fs);
