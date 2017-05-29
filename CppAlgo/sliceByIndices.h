#pragma once
#include <Eigen/Dense>
#include <type_traits>
#include <vector>
#include "types.h"

/**
 * Slice a vector by a given index array
 *@tparam Derived the input vector type
 *@tparam DerivedIndex the index vector type
 *@param v the input vector
 *@param indices the index vector 
 *@param indexBase whether the given indices are one-based or zero-based, default to be zero-based
 *@return the index-sliced subvector
 *@note If we slice a vector [-1, 2, -3, 4, 5] by a index array (one-based) [1, 3, 4], we get [-1, -3, 4]
*/
template<typename Derived, typename DerivedIndex>
typename Derived::PlainObject sliceByIndices(const Eigen::MatrixBase<Derived>& v, const Eigen::MatrixBase<DerivedIndex>& indices, IndexBase indexBase = IndexBase::Zero)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedIndex);
	static_assert(std::is_integral<typename DerivedIndex::Scalar>::value, "The type of indices must be integers.");
	auto size = indices.size();
	typename Derived::PlainObject r(size);
	int diff = indexBase == IndexBase::One ? 1 : 0;
	for (int i = 0; i < size; i++)
		r(i) = v(indices(i) - diff);
	return r;
}



/**
 * Index a vector by logicals.
 *@tparam Derived type of the source
 *@tparam DerivedLogical type of the logical (its scalar type should be integral)
 *@tparam Result type of the result
 *@param source input vector/array to be indexed
 *@param logical input logical vector/array for indexing
 *@return the indexed vector/array
 **@note this is used to mimic the MATLAB logical array indexing infrastructure, see https://www.mathworks.com/help/matlab/math/matrix-indexing.html?refresh=true
*/
template<typename Derived, typename DerivedLogical, typename Result = typename Derived::PlainObject>
Result indexByLogical(const Eigen::DenseBase<Derived>& source, const Eigen::DenseBase<DerivedLogical>& logical)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedLogical);
	static_assert(std::is_integral<typename DerivedLogical::Scalar>::value, "The logical must be of integral types.");

	Result r(source.size());
	Eigen::Index j = 0;
	for (Eigen::Index i = 0; i < logical.size(); i++)
	{
		if (logical(i))
		{
			r(j) = source(i);
			++j;
		}
	}
	r.conservativeResize(j);
	return r;
}
