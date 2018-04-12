#pragma once
#include<Eigen/Dense>
#include "types.h"

//@brief A structure type for picos
struct PicosElement
{
	int pm;
	Eigen::TFloat f0;
	Eigen::TRowVectorX a;
	Eigen::TRowVectorX p;
	Eigen::TRowVectorX e;
	Eigen::TRowVectorX mfcc;
	Eigen::TFloat alfa;
};