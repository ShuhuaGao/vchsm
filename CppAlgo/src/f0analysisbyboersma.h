#pragma once
#include <Eigen/Dense>
#include "types.h"

Eigen::TRowVectorX f0analysisbyboersma(const Eigen::Ref<const Eigen::TRowVectorX>& x, Eigen::TFloat fs, const Eigen::Ref<const Eigen::RowVectorXi>& pms,
	Eigen::TFloat f0min = 60, Eigen::TFloat f0max = 500);
