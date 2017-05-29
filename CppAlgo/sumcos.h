#pragma once
#include<Eigen/Dense>
#include "types.h"

enum class Modo
{
	sin,
	cos
};

Eigen::TRowVectorX sumcos(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::Ref<const Eigen::TRowVectorX> pp, Modo modo);