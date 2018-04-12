#pragma once
#include "types.h"

template<typename FloatType>
struct Tolerance
{
	using type = FloatType;
};

template<>
struct Tolerance<float>
{
	using type = float;
	constexpr static float absoluteTol = 1e-5F;
	constexpr static float relativeTol = 1e-6F;
};

template<>
struct Tolerance<double>
{
	using type = double;
	constexpr static double absoluteTol = 1e-12;
	constexpr static double relativeTol = 1e-12;
};

/**
 * Check whether two floating numbers are approximately equal
 *tparam FloatType floating type for the comparison, float or double (corresponding to different tolerances)
 *@param a the first number
 *@param b the second number
 *@return true - equal, false - not equal
*/
template<typename FloatType = Eigen::TFloat, typename T1, typename T2>
bool almostEqual(T1 a, T2 b, FloatType absoluteTol = Tolerance<FloatType>::absoluteTol, 
	FloatType relativeTol = Tolerance<FloatType>::relativeTol)
{
	// first check the absolute tolerance
	auto absDiff = std::abs(a - b);
	if (absDiff <= absoluteTol)
		return true;
	// then check the relative tolerance
	auto abs_a = std::abs(a);
	auto abs_b = std::abs(b);
	if (absDiff <= (abs_a > abs_b ? abs_a : abs_b) * relativeTol)
		return true;
	return false;
}
