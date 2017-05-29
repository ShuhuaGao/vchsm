#include "stochasticanalysis.h"
#include "seq.h"
#include "constants.h"
#include "filter.h"
#include "bymcaquat.h"
#include "concat.h"
#include "bylevdurb.h"
#include "almostEqual.h"
#include <limits>

// Clang cannot link Tolerance<double>, which is quite OK in visual studio.
// Just define it as a macro here
#define Tolerance_double 1.0e-12


void stochasticanalysis(Eigen::Ref<const Eigen::TRowVectorX> x, Eigen::TFloat fs, int N, PicosStructArray& picos, int ordenLPC)
{
	int Npm = static_cast<int>(picos.size());
	Eigen::TRowVectorX yd(x.size());
	yd.setZero();

	auto left = picos[0].pm - N - 1; // segment starts at left and length is N
	auto size = picos[0].a.size();
	Eigen::TRowVectorX seq1 = seq(0, N - 1) / Eigen::TFloat(N); // (0:N-1)/N
	Eigen::TRowVectorX seq2 = (2 * pi / fs * picos[0].f0 ) *seq(-N, -1); // 2*pi*picos(1).f0*(-N:-1)/fs
	for (int j = 1; j <= size; j++)
	{
		yd.segment(left, N) += ((picos[0].a(j - 1) * seq1).array() * (j * seq2.array() + picos[0].p(j - 1)).cos()).matrix();
	}

	Eigen::TRowVectorX seq3 = seq(0, N - 1); // 0: N-1
	Eigen::TRowVectorX seq4 = seq(1, N).reverse(); // N:-1:1
	auto seq5 = seq(-N, -1); //-N: -1
	for (int k = 1; k <= Npm - 1; k++)
	{
		auto minSize = std::min(picos[k - 1].a.size(), picos[k].a.size());
		for (int j = 1; j <= minSize; j++)
		{
			auto a1 = picos[k - 1].a(j - 1);
			auto a2 = picos[k].a(j - 1);
			Eigen::TFloat c0, c1, c2, c3;
			std::tie(c0, c1, c2, c3) = bymcaquat(picos[k - 1].p(j - 1), picos[k].p(j - 1), j * picos[k - 1].f0, j * picos[k].f0, N, fs);
			yd.segment(picos[k - 1].pm - 1, N).array() += (a1 + ((a2 - a1) / N) * seq3.array()) * (c0 + seq3.array() * (c1 + seq3.array() * (c2 + c3 * seq3.array()))).cos();
		}


		for (int j = (int)picos[k].a.size() + 1; j <= picos[k-1].a.size(); j++)
		{
			yd.segment(picos[k - 1].pm - 1, N).array() += (picos[k - 1].a(j - 1) / N * seq4.array()) 
				* ((2 * pi*j*picos[k - 1].f0 / fs) * seq3.array() + picos[k - 1].p(j - 1)).cos();
		}

		for (int j = (int)picos[k - 1].a.size() + 1; j <= picos[k].a.size(); j++)
		{
			yd.segment(picos[k].pm - N - 1, N).array() += (picos[k].a(j - 1) / N * seq3).array() * (2 * pi*j*picos[k].f0 / fs * seq5.array() + picos[k].p(j - 1)).cos();
		}

	}

	// to avoid copying ye, we put the following zeros here 
	Eigen::TRowVectorX ye = concat(x - yd, Eigen::TRowVectorX::Zero(24));
	// ye=filter(fir1(48,500/(0.5*fs),'high'),1,ye); MATLAB
	// It is found that fs is fixed to be 16000Hz. Therefore fir1(48,500/(0.5*fs),'high') returns a fixed value, which exempts us from implementation of fir1.
	filter(ye);
	ye.tail(24).setZero();
	auto ye_ = ye.tail(ye.size() - 24);
	Eigen::TRowVectorX hnnN = std::sqrt(8.0 / 3.0) * (0.5 - 0.5*(2 * pi*seq1).array().cos());

	// better to allocate the fixed-length matrix outside for-loop
	Eigen::TRowVectorX trama(N);
	Eigen::TRowVectorX R(1 + ordenLPC);
	for (int k = 1; k <= Npm; k++)
	{
		trama = ye_.segment(picos[k - 1].pm - static_cast<int>(std::floor(N / 2.0)) - 1, N).cwiseProduct(hnnN);
		if (trama.isZero(Tolerance_double))
		{
			picos[k - 1].e = Eigen::TRowVectorX::Zero(ordenLPC + 1);
			picos[k - 1].e(0) = std::numeric_limits<Eigen::TFloat>::infinity();
		}
		else
		{
			R.setZero();
			for (int j = 0; j <= ordenLPC; j++)
			{
				R(j) = trama.tail(N - j).dot(trama.head(N - j));
			}
			picos[k - 1].e = std::sqrt(N) * bylevdurb(R);
		}
			
	}

}
