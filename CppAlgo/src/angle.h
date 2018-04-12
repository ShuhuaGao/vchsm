#pragma once
#include <Eigen/Dense>
#include "types.h"

/**
  * Compute the phase angle for each element of a complex matrix.
  *@param m input matrix
  *@return angle matrix storing the angle of each element in m
  *@note Implementation of MATLAB angle function. see https://www.mathworks.com/help/matlab/ref/angle.html
  *@remarks Mathematically, angle(z) = imag(log(z)) = atan2(imag(z),real(z)).
*/
template<typename Derived>
auto angle(const Eigen::MatrixBase<Derived>& m)
{
	auto findAngle = [](const std::complex<Eigen::TFloat>& c) {return std::atan2(c.imag(), c.real());};
	return m.unaryExpr(findAngle).real(); // apply unitary operation to matrix m coefficient-wise.
}


/**
* Compute the phase angle for each element of a complex number
*@param c input scalar complex number
*@return angle 
*@note Implementation of MATLAB angle function. see https://www.mathworks.com/help/matlab/ref/angle.html
*@remarks Mathematically, angle(z) = imag(log(z)) = atan2(imag(z),real(z)).
*/
template<typename T>
auto angle(const std::complex<T>& c)
{
	return std::atan2(c.imag(), c.real());
}
