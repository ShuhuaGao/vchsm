#pragma once
#include <cmath>

/**
 * Round up to next higher power of 2 (return x if it's already a power)
 *@param x  a integer
 *@return a value y = 2^p such that y is the smallest powers of two that is satisfies 2^p >= |x|
 *@remark This is a similar implementation of MATLAB nextpow2 function. However, note that the MATLAB one will return p instead of y.
 *@note This code is copied from http://stackoverflow.com/questions/364985/algorithm-for-finding-the-smallest-power-of-two-thats-greater-or-equal-to-a-giv.
*/
int pow2Roundup(int x)
{
	x = std::abs(x);
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}

