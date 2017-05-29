#include "HSMsynthesize.h"
#include "seq.h"
#include "concat.h"
#include "constants.h"
#include <algorithm>
#include <random>

// filter(1,picos(k).e,randn(1,n3-n1))
// see https://www.mathworks.com/help/matlab/ref/filter.html
// b = 1, a = picos(k).e, x = randn(1,n3-n1)
// a(1)y(n) = x(n) - a(2)y(n-1) -.. - a(na+1)y(n-na)
// if we first normalize a (i.e., pe), then y(n) = 1/a(1)*x(n) - a(2)/a(1)*y(n-1) - ...- a(na+1)/a(1)*y(n-na)
Eigen::TRowVectorX filter(const Eigen::TRowVectorX& a, const Eigen::TRowVectorX& x)
{
	auto na = a.size() - 1;
	Eigen::TRowVectorX an = a / a(0); // normalize
	Eigen::TRowVectorX y(x.size());
	for (Eigen::Index i = 0; i < na; i++)
	{
		y(i) = x(i) / a(0) - an.segment(1, i).dot(y.head(i).reverse());
	}
	for (Eigen::Index  i = na; i < y.size(); i++)
	{
		y(i) = x(i) / a(0) - an.segment(1, na).dot(y.segment(i - na, na).reverse());
	}
	return y;
}


Eigen::TRowVectorX synth(int L, int fs, const PicosStructArray & audioFeature)
{
	assert(fs == 16000);
	const auto& picos = audioFeature;
	Eigen::TRowVectorX y(L);
	y.setZero();

	auto Npm = picos.size();

	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<Eigen::TFloat> nd; // standard normal distribution.

	for (std::size_t k = 1; k <= Npm; k++)
	{
		auto n2 = picos[k - 1].pm;
		int n1 = 0, n3 = 0;
		if (k == 1)
		{
			n3 = picos[k].pm;
			n1 = std::max(1, n2 - (n3 - n2));
		}
		else if (k == Npm)
		{
			n1 = picos[k - 2].pm;
			n3 = std::min(L, n2 + n2 - n1);
		}
		else
		{
			n1 = picos[k - 2].pm;
			n3 = picos[k].pm;
		}
		auto N12 = n2 - n1;
		auto N23 = n3 - n2;
		auto seq1 = seq(0, 1.0 / N12, (N12 - 1) / (Eigen::TFloat)N12);
		auto seq2 = seq(1.0 / N23, 1.0 / N23, 1).reverse();
		auto win = concat<Eigen::TRowVectorX>(seq1, seq2);
		Eigen::TRowVectorX trama(n3 - n1);
		trama.setZero();
		for (Eigen::Index j = 1; j <= picos[k - 1].a.size(); j++)
		{
			trama += picos[k - 1].a(j - 1) * (2 * pi * seq(-N12, N23 - 1).array() * (j * picos[k - 1].f0) / fs + picos[k - 1].p(j - 1) + j * picos[k - 1].alfa).cos().matrix();
		}
		// randn(1,n3-n1), generate the normally distributed random vector
		
		Eigen::TRowVectorX x(n3 - n1);
		for (Eigen::Index i = 0; i < x.size(); ++i)
		{
			x(i) = nd(gen);
		}

		trama += filter(picos[k - 1].e, x);
		y.segment(n1 - 1, n3 - n1) += trama.cwiseProduct(win);
	}
	return y;
}


Eigen::TRowVectorX HSMsynthesize(const PicosStructArray & audioFeature, int L)
{
	constexpr int fs = 16000;
	return synth(L, fs, audioFeature);
}
