#pragma once
#include <Eigen/Dense>
#include <assert.h>
#include "types.h"

/**
 * Concatenate two compatible vectors or matrices
 *@tparam ResultType type of the constructed concatenation vector or matrix, default to be Eigen::TRowVectorX
 *@tparam Derived1 type of the 1st vector or matrix to be concatenated
 *@tparam Derived2 type of the 2nd vector or matrix to be concatenated
 *@param X first vector or matrix
 *@param Y second vector or matrix
 *@param direction the direction for concatenation, default to be horizontal
 *@return the concatenated one 
 *@remarks For example, if X and Y are two row vectors, then concat(X, Y, Horizontal) will behave like [X, Y] in MATLAB
 *@note Derived1 and Derived2 must have the same Scalar type, such as Eigen::TFloat, int or float, .etc.
*/
template<typename ResultType = Eigen::TRowVectorX, typename Derived1, typename Derived2>
auto concat(const Eigen::MatrixBase<Derived1>& X, const Eigen::MatrixBase<Derived2>& Y,
	Eigen::DirectionType direction = Eigen::Horizontal)
{
	assert(direction == Eigen::Horizontal && X.rows() == Y.rows()
		|| direction == Eigen::Vertical && X.cols() == Y.cols());
	if (direction == Eigen::Horizontal)
	{
		return (ResultType(X.rows(), X.cols() + Y.cols()) << X, Y).finished();
	}
	else
	{
		return (ResultType( X.rows() + Y.rows(), X.cols()) << X, Y).finished();
	}
	
}

/**
 * Concatenate a vector and a scalar 
*/
template<typename ResultType = Eigen::TRowVectorX, typename Derived>
auto concat(const Eigen::MatrixBase<Derived>& X, const typename Derived::Scalar& y, Eigen::DirectionType direction = Eigen::Horizontal)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(ResultType);
	assert(direction == Eigen::Horizontal && X.rows() == 1 ||
		direction == Eigen::Vertical && X.cols() == 1);
	if (direction == Eigen::Horizontal)
	{
		return (ResultType(X.cols() + 1) << X, y).finished();
	}
	else
	{
		return (ResultType(X.rows() + 1) << X, y).finished();
	}
}


/**
* Concatenate a scalar x and a vector Y to form [x, Y]
*/
template<typename ResultType = Eigen::TRowVectorX, typename Derived>
auto concat(const typename Derived::Scalar& x, const Eigen::MatrixBase<Derived>& Y, Eigen::DirectionType direction = Eigen::Horizontal)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(ResultType);
	assert(direction == Eigen::Horizontal && Y.rows() == 1 ||
		direction == Eigen::Vertical && Y.cols() == 1);
	if (direction == Eigen::Horizontal)
	{
		return (ResultType(Y.cols() + 1) << x, Y).finished();
	}
	else
	{
		return (ResultType(Y.rows() + 1) << x, Y).finished();
	}
}




