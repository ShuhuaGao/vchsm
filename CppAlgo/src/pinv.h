#pragma once

#include <Eigen/Dense>

/**
 * Moore-Penrose pseudoinverse of matrix. Implementation of MATLAB inv function.
 *@tparam Derived the input matrix type
 *@param mat input matrix
 *@param tolerance every singluar value of mat that is less than tolerance will be treated as zero.
 *@return the pseudoinverse
 *@note This code is adapted from https://gist.github.com/javidcf/25066cf85e71105d57b6.
*/
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::ColsAtCompileTime, Derived::RowsAtCompileTime>
pinv(const Eigen::MatrixBase<Derived> &mat, typename Derived::Scalar tolerance = Eigen::NumTraits<typename Derived::Scalar>::dummy_precision()) 
{
	typedef typename Derived::Scalar Scalar;
	auto svd = mat.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
	const auto &singularValues = svd.singularValues();
	Eigen::Matrix<Scalar, Derived::ColsAtCompileTime, Derived::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
	singularValuesInv.setZero();
	for (int i = 0; i < singularValues.size(); ++i) {
		if (singularValues(i) > tolerance)
		{
			singularValuesInv(i, i) = Scalar{ 1 } / singularValues(i);
		}
		else
		{
			singularValuesInv(i, i) = Scalar{ 0 };
		}
	}
	return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}