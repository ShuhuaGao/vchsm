#pragma once
#include <Eigen/Dense>
#include "types.h"

Eigen::TFloat linphaseterm(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::Ref<const Eigen::TRowVectorX> pp, int fmaxopt = 1000);
