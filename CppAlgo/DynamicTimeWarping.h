#pragma once
#include "types.h"
#include <utility>
#include <vector>

std::pair<std::vector<int>, std::vector<int>> DynamicTimeWarping(Eigen::Ref<const Eigen::TMatrixX> A, Eigen::Ref<const Eigen::TMatrixX> B);
