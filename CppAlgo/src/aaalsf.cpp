#include "aaalsf.h"
#include "seq.h"
#include "constants.h"
#include "angle.h"
#include "concat.h"
#include "roots.h"
#include <algorithm>

Eigen::TVectorX aaalsf(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::TFloat f0, int p)
{
	const int fmax = 5000;
	int Nk = (int)aa.size();
	auto ff = (seq(1, Nk) * f0).eval();
	auto PP = aa.array().square().eval();
	auto fs = 2 * fmax;
	auto R = Eigen::TRowVectorX::Zero(p + 1).eval();
	for (int j = 1; j <= p + 1; j++)
	{
		R(j - 1) = (1.0 / Nk) * (PP * (2 * pi*(j - 1) / fs * ff).array().cos()).sum();
	}
	Eigen::TRowVectorX ai(p);
	ai.setZero();
	auto e = R(0);
	for (int j = 1; j <= p; j++)
	{
		if (j == 1)
		{
			auto k = R(1) / R(0);
			ai(0) = k;
			e = (1 - k*k) * e;
		}
		else
		{
			auto k = R(j);
			k = k - ai.head(j - 1).dot(R.segment(1, j - 1).reverse());
			k = k / e;
			ai(j - 1) = k;
			Eigen::TRowVectorX aux = ai.head(j - 1) - k * ai.head(j - 1).reverse();
			ai.head(j - 1) = aux;
			e = (1 - k* k) * e;
		}
	}
	Eigen::TRowVectorX az1(p + 2); // length of az1 is 2 + length(ai)
	//auto temp = std::sqrt(std::abs(R(0) - ai.dot(R.tail(p)))); // sqrt(abs( R(1)-sum(ai.*R(2:p+1)) ))
	az1(0) = 1;
	az1.segment(1, p) = -ai;
	az1(p + 1) = 0;
	auto az2 = az1.reverse();
	auto lsf = angle(concat<Eigen::TRowVectorXc>(roots(az1 + az2).transpose(), roots(az1 - az2).transpose())).eval();
	int lsfSize = (int)lsf.size();
	std::sort(lsf.data(), lsf.data() + lsfSize);

	return lsf.segment(lsfSize - p - 1, p).transpose();
}
