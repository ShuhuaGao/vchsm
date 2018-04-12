#pragma once

#include <utility>
#include <Eigen/Dense>
#include "types.h"

using Eigen::TMatrixX;
using Eigen::TRowVectorX;

std::pair<TMatrixX, TMatrixX> polytseb(int K);
