#pragma once
#include "types.h"
#include "almostEqual.h"

/**
 * Find polynomial roots. An implementation of roots function of MATALB.
 *@tparam Derived the polynomial coefficient list type
 *@param c the polynomial coefficient list
 *@return the complex roots in a column vector
 *@note Use "edit roots" to scan the source code in MATALB.
*/
template<typename Derived>
auto roots(const Eigen::MatrixBase<Derived>& c)
{
	using Scalar = typename Derived::Scalar;
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	auto n = (int)c.size();

	// Strip leading zeros and throw away.  
	// Strip trailing zeros, but remember them as roots at zero.
	int i = 0;
	for (i = 0; i < n; ++i)
	{
		if (!almostEqual(c(i), 0))
			break;
	}
	int j = 0;
	for (j = n - 1; j >= 0; j--)
	{
		if (!almostEqual(c(j), 0))
			break;
	}
	// now c(i:j) is not zero (i and j included)
	int nr = n - j - 1; // number of trailing zeros
	using ComplexVectorType = Eigen::Matrix<std::complex<Scalar>, -1, 1>;
	auto r = ComplexVectorType::Zero(nr, 1).eval();
	// Prevent relatively small leading coefficients from introducing Inf
	// by removing them
	int k = i + 1;
	while (k <= j && !(c.segment(k, j - k + 1) / c(k - 1)).allFinite())
	{
		k++;
	}
	auto d = c.segment(k, j - k + 1) / c(k - 1);

	// Polynomial roots via a companion matrix
	int nn = j - k + 2;
	if (nn > 1)
	{
		auto a = Eigen::Matrix<Scalar, -1, -1>::Zero(nn - 1, nn - 1).eval();
		for (int i = 1; i < a.rows(); i++)
		{
			a(i, i - 1) = 1;
		}
		a.row(0) = -d;
		auto eig = Eigen::EigenSolver<Eigen::Matrix<Scalar, -1, -1>>(a, false).eigenvalues();
		ComplexVectorType cv(r.size() + eig.size());
		cv << r, eig;
		return cv;
	}
	return r;	
}
