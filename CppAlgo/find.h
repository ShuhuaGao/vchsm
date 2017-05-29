#pragma once
#include "types.h"

/**
 * Find linear indices of nonzero elements in a vector or matrix.
 *@tparam Derived type of the input matrix/vector
 *@param x the input matrix/vector
 *@param indexBase the indexing base of the returned indices, default to be 0
 *@return a column vector containing the indices of nonzero elements in x. 
 *@note If x is a matrix,  this function returns a column vector of the linear indices of the result.
 *@remark This is an implementation of find function in MATLAB. See https://www.mathworks.com/help/matlab/ref/find.html.
*/
template<typename Derived>
Eigen::Matrix<Eigen::Index, -1, 1> find(const Eigen::MatrixBase<Derived>& x, IndexBase indexBase = IndexBase::Zero)
{
	Eigen::Matrix<Eigen::Index, -1, 1> r(x.size());
	Eigen::Index j = 0;
	int diff = indexBase == IndexBase::Zero ? 0 : 1;
	for (Eigen::Index i = 0; i < x.size(); i++)
	{
		// Eigen also supports linear indexing (column-major as default)
		if (x(i) != 0)
		{
			r(j) = i + diff;
			j++;
		}
	}
	r.conservativeResize(j); // this will succeed even if j is 0 --> an empty one
	return r;
}
