#pragma once
#include <tuple>
#include <Eigen/Dense>
#include "types.h"

std::tuple<Eigen::TFloat, Eigen::TFloat, Eigen::TFloat, Eigen::TFloat> bymcaquat(Eigen::TFloat p1, Eigen::TFloat p2, Eigen::TFloat f1, Eigen::TFloat f2, Eigen::TFloat N, Eigen::TFloat fs);
