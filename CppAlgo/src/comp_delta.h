#pragma once
#include <Eigen/Dense>
#include "types.h"

/**
 * Implementation of comp_delta in comp_delta.m.
*/
Eigen::TMatrixX comp_delta(const Eigen::Ref<const Eigen::TMatrixX>& static_coef, int DELTAWINDOW);

