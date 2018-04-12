#include "sumcos.h"
#include "polytseb.h"
#include "conv.h"
#include "roots.h"
#include "polyeval.h"
#include "almostEqual.h"
#include <algorithm>
#include <vector>

Eigen::TRowVectorX findRealAndMagLessThanOne(const Eigen::TVectorXc& x);

Eigen::TRowVectorX sumcos(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::Ref<const Eigen::TRowVectorX> pp, Modo modo)
{
	auto K = (int)aa.size();
	auto AcP = aa.cwiseProduct(pp.array().cos().matrix()); // cannot mix matrix and array, convert to matrix back
	auto AsP = aa.cwiseProduct(pp.array().sin().matrix());
	Eigen::TMatrixX T, U;
	std::tie(T, U) = polytseb(K);

	Eigen::TRowVectorX Tx, Ux;
	if (modo == Modo::cos)
	{
		Tx.noalias() = AcP * T;
		Ux.noalias() = -AsP * U;
	}
	else
	{
		Ux.noalias() = AcP * U;
		Tx.noalias() = AsP * T;
	}

	auto UxConv = conv(Ux, Ux);
	auto TxConv = conv(Tx, Tx);
	Eigen::TRowVectorX Px = conv((Eigen::TRowVectorX(3) << -1, 0, 1).finished(), UxConv)
		- (Eigen::TRowVectorX(2 + TxConv.size()) << 0, 0, TxConv).finished();

	auto x = findRealAndMagLessThanOne(roots(Px));

	Eigen::TFloat term1, term2, xk;
	for (int k = 1; k <= x.size(); k++)
	{
		xk = x(k - 1); // index from 0
		term1 = std::sqrt(1 - xk*xk) * polyval(Ux, xk);
		term2 = polyval(Tx, xk);
		if (std::abs(term2 + term1) < std::abs(term2 - term1))
		{
			x(k - 1) = std::acos(xk);
		}
		else
		{
			x(k - 1) = -std::acos(xk);
		}
	}

	return x;
}

//// find elements in x that are real and magnitude <= 1
Eigen::TRowVectorX findRealAndMagLessThanOne(const Eigen::TVectorXc& x)
{
	std::vector<Eigen::TFloat> v;
	v.reserve(x.size());
	for (int i = 0; i < x.size(); ++i)
	{
		if (almostEqual(x(i).imag(), 0) && std::abs(x(i)) <= 1)
		{
			v.push_back(x(i).real());
		}
	}
	Eigen::TRowVectorX r(v.size());
	std::copy(v.cbegin(), v.cend(), r.data());
	return r;
}