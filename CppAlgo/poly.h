#pragma once
#include "types.h"
#include <unsupported/Eigen/Polynomials>

/**
 * Find the Polynomial with specified roots.
 *@tparam Derived vector type of the roots
 *@param r the roots 
 *@return a row vector containing the polynomial coefficients in dscending order w.r.t to the degree
 *@remark This is an implementation of MATLAB poly function, see https://www.mathworks.com/help/matlab/ref/poly.html
 *@remark If the input roots are complex, the polynomial coefficients may also be complex. If in fact the coefficients should be real, use 
 * real() method to get the real parts.
*/
template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, 1, -1> poly(const Eigen::MatrixBase<Derived>& r)
{
	Eigen::Matrix<typename Derived::Scalar, 1, -1> polynomial;
	Eigen::roots_to_monicPolynomial(r, polynomial);
	return polynomial.reverse();
}