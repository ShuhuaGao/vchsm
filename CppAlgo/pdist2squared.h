#pragma once
#include <Eigen/Dense>
#include "types.h"

/**
 * Pairwise distance between two sets of observations. Implementation of MATLAB pdist2 function.
 *@tparam Derived A matrix type
 *@param X first matrix
 *@param Y second matrix
 *@return the pairwise distance (L2 norm) matrix
 *@details D = pdist2(X,Y) returns a matrix D containing the Euclidean distances between each pair of observations in the mx-by-n data matrix X and my-by-n data matrix Y. 
 * Rows of X and Y correspond to observations, columns correspond to variables. D is an mx-by-my matrix, with the (i,j) entry equal to distance between observation i in X and observation j in Y. The (i,j) entry will be NaN if observation i in X or observation j in Y contain NaNs.
 *@note see https://www.mathworks.com/help/stats/pdist2.html
 */
template<typename Derived>
typename Derived::PlainObject pdist2(const Eigen::MatrixBase<Derived>& X, const Eigen::MatrixBase<Derived>& Y)
{
	typename Derived::PlainObject D(X.rows(), Y.rows());
	for (int i = 0; i < X.rows(); i++)
	{
		D.row(i) = (Y.rowwise() - X.row(i)).rowwise().norm();
	}
	return D;
}


/**
* Squared pairwise distance between two sets of observations.
*@see pdist2
*/
template<typename Derived>
typename Derived::PlainObject pdist2Squared(const Eigen::MatrixBase<Derived>& X, const Eigen::MatrixBase<Derived>& Y)
{
	typename Derived::PlainObject D(X.rows(), Y.rows());
	for (int i = 0; i < X.rows(); i++)
	{
		D.row(i) = (Y.rowwise() - X.row(i)).rowwise().squaredNorm();
	}
	return D;
}


