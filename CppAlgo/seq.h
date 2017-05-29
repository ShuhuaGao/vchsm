#pragma once
#include <Eigen/Dense>

/**
  *Generate a sequence (vector) like the MATALB colon operation
  *@tparam Derived a row or column vector type, default to be Eigen::TRowVectorX
  *@param start number, included
  *@param step size between every adjacent numbers, must be positive
  *@param end number, may not be included, depends on start and step
  *@return an Eigen::CwiseNullaryOp expression, which can used like a vector
  *@note seq(1, 1, 5) is like 1:1:5 (or 1:5 for short) which gives [1, 2, 3, 4, 5]
  *@note seq(start, step, end) is like start:step:end in MATLAB
  *@remarks see https://www.mathworks.com/help/matlab/ref/colon.html
  *
*/
template<typename Derived = Eigen::TRowVectorX>
auto seq(typename Derived::Scalar start, typename Derived::Scalar step, typename Derived::Scalar end)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	assert(start <= end); // must ensure start <= end
	int m = (int)std::trunc((end - start) / step); // std::trunc <--> fix() in MATLAB
	return Eigen::DenseBase<Derived>::LinSpaced( m + 1, start, start + m*step);
}

/**
*Generate a sequence (vector) like the MATALB colon operation with a default step size = 1
*@tparam Derived a row or column vector type, default to be Eigen::TRowVectorX
*@param start number, included
*@param end number, may not be included, depends on start and step = 1
*@return an Eigen::CwiseNullaryOp expression, which can used like a vector
*@note seq(1,5) is like 1:5 (or 1:1:5) which gives [1, 2, 3, 4, 5]
*@note seq(start, end) is like start:end in MATLAB
*@remarks see https://www.mathworks.com/help/matlab/ref/colon.html
*@sa seq(Derived::Scalar start, Derived::Scalar step, Derived::Scalar end)
*/
template<typename Derived = Eigen::TRowVectorX>
auto seq(typename Derived::Scalar start, typename Derived::Scalar end)
{
	return seq<Derived>(start, 1, end);
}