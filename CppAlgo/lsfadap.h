#pragma once
#include "types.h"

/**
 * Implementation of lsfadap.m.
 *@note the 2nd input argument factor in lsfadap.m is never used, which is therefore neglected here.
*/
Eigen::TRowVectorX lsfadap(Eigen::Ref<const Eigen::TVectorX> lsf);