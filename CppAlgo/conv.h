#pragma once
#include <Eigen/Dense>
#include <algorithm>
using Eigen::Index;

/** Vector convolution
  *
  * \param u one vector (row or column)
  * \param v another vecotr (row or column)
  * \return the convolution of vectors u and v (row or column, same as u)
  * \note This is the implementation of MATLAB conv function. Check https://www.mathworks.com/help/matlab/ref/conv.html#bucr92l-2.
  * \remark Use template here since we want to support both row and column vector inputs.
  */
template<typename TVector1, typename TVector2>
TVector1 conv(const TVector1 & u, const TVector2 & v)
{
	Index m = u.size();
	Index n = v.size();
	// the convolution result is a vector of length m + n - 1
	TVector1 w(m + n - 1);
	// compute each element of w
	for (Index k = 1; k <= m + n - 1; k++)
	{
		Index jMin = std::max<Index>(1, k + 1 - n);
		Index jMax = std::min<Index>(k, m);
		// Note: Eigen indexes from 0
		w(k - 1) = u.segment(jMin - 1, jMax - jMin + 1).dot(v.segment(k - jMax, jMax - jMin + 1).reverse());
	}
	return w;
}


