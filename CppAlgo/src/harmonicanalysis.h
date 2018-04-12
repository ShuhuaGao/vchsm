#pragma once

#include <Eigen/Dense>
#include <vector>
#include "PicosElement.h"
#include "types.h"
#include "StructDefinitions.h"


PicosStructArray harmonicanalysis(Eigen::Ref<const Eigen::TRowVectorX>  x, Eigen::TFloat fs, Eigen::Ref<const Eigen::RowVectorXi> pms, Eigen::Ref<const Eigen::TRowVectorX> f0s, Eigen::TFloat fmax);
