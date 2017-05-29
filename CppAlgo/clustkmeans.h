#pragma once
#include "types.h"
#include <utility>

std::pair<Eigen::RowVectorXi, Eigen::RowVectorXi> clustkmeans(Eigen::Ref<const Eigen::TMatrixX> X, int n);