#pragma once
#include <unsupported/Eigen/Polynomials>

/**
  * Find the value of the given polynomial at x
  * @param p the polynomial coefficient vector in descending power order
  * @param x a scalar value
  * @tparam VectorP a vector type
  * @tparam T a scalar type
  * @return value of the polynomial at x
  * @remark Implementation of MATLAB polyval, but only for a scalar input.
  * @remark polyval([1, 2, 3], 10.1) finds x^2 + 2*x + 3 at x = 10.1
  */
template<typename Derived, typename T>
T polyval(const Eigen::MatrixBase<Derived>& p, T x)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	return Eigen::poly_eval(p.reverse(), x);
}