#pragma once
#include <Eigen/Dense>
#include "seq.h"
#include "constants.h"

/**
 * Hamming window
 *@param L length of the Hamming window
 *@return a L-point sysmmetric Hamming window in the column vector
 *@ramark implementation of the hamming function of MATLAB, see https://www.mathworks.com/help/signal/ref/hamming.html
*/
Eigen::TVectorX hamming(int L)
{
	int N = L - 1;
	return 0.54 - 0.46 * (2 * pi / Eigen::TFloat(N)  * seq<Eigen::TVectorX>(0, N)).array().cos();
}