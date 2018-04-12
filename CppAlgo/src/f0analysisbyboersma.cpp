#include "f0analysisbyboersma.h"
#include "seq.h"
#include "constants.h"
#include "FFT.h"
#include "concat.h"
#include "almostEqual.h"
#include "sliceByIndices.h"
#include <vector>
#include <limits>

/**
 *@struct DatElement
 *Structure type for dat
 *@note MATLAB: dat(1:length(pms))=struct('r',[],'f',[],'ac',[],'ant',[]);
*/
struct DatElement
{
	Eigen::TRowVectorX r, f, ac;
	Eigen::RowVectorXi ant;
};

Eigen::TRowVectorX f0analysisbyboersma(const Eigen::Ref<const Eigen::TRowVectorX>& x, Eigen::TFloat fs, const Eigen::Ref<const Eigen::RowVectorXi>& pms,
	Eigen::TFloat f0min, Eigen::TFloat f0max )
{
	int Ncandidates = 6;
	int  lagmin = static_cast<int>(std::ceil(fs / f0max));
	int lagmax = static_cast<int>(std::floor(fs / f0min));
	Eigen::TFloat fact = 0.01 * fs / (pms(1) - pms(0));
	Eigen::TFloat voith = 0.45;
	Eigen::TFloat silth = 0.03;
	Eigen::TFloat octcost = 0.01;
	Eigen::TFloat vuvcost = 0.14*fact;
	Eigen::TFloat uvvcost = 0.14*fact;
	Eigen::TFloat uvuvcost = 0.0;
	Eigen::TFloat octjump = 0.35*fact;
	int L = static_cast<int>(std::ceil(3.0*fs / f0min));
	int L2 = 0;
	if (L % 2 == 0)
	{
		L2 = L / 2;
		L++;
	}
	else
	{
		L2 = (L - 1) / 2;
	}

	int Lz = static_cast<int>(std::ceil(1.5 * L));
	int Lp2 = 2;
	while (Lp2 < Lz)
	{
		Lp2 *= 2;
	}

	int Lxlmax = static_cast<int>(std::ceil(0.5 * fs / f0min));
	int Ldc = 2 * Lxlmax;
	Eigen::TRowVectorX w = Eigen::TRowVectorX::Zero(Lp2);
	w.head(L) = 0.5 * (1 - (2 * pi * seq(0, L - 1).array() / L).cos());

	//TODO: change the FFT to iOS vDisp implementation
	auto tt = fft(w).array().abs2().matrix().cast<std::complex<Eigen::TFloat>>().eval(); // must be complex type though the number is actually real
	auto rw_temp = ifft(tt).real().eval(); // cannot use expression template here since it will involve the temporary returned from the function

	auto rw = (1 / rw_temp(0)) * rw_temp.head(L2);
	Eigen::TFloat xgmax = x.cwiseAbs().maxCoeff();
	std::vector<DatElement> dat(pms.size());
	Eigen::TRowVectorX rx(L2);
	Eigen::TRowVectorX rk(Ncandidates - 1);
	Eigen::TRowVectorX fk(Ncandidates - 1);

	for (Eigen::Index k = 1; k <= pms.size(); k++)
	{
		Eigen::TRowVectorX trama;
		if (pms(k - 1) - L2 < 1)
		{
			auto tramaTemp = x.head(pms(k - 1) + L2);
			trama = concat(Eigen::TRowVectorX::Zero(L - tramaTemp.size()), tramaTemp);
		}
		else if (pms(k - 1) + L2 > x.size())
		{
			int length = (int)x.size() - (pms(k - 1) - L2) + 1;
			auto tramaTemp = x.tail(length);
			trama = concat(tramaTemp, Eigen::TRowVectorX::Zero(L - tramaTemp.size()));
		}
		else
		{
			trama = x.segment(pms(k - 1) - L2 - 1, 2 * L2 + 1);
		}

		trama.array() -= trama.segment(L2 - Ldc, 2 * Ldc + 1).sum() / (2 * Ldc + 1.0);
		Eigen::TFloat xlmax = trama.segment(L2 - Lxlmax, 2 * Lxlmax + 1).cwiseAbs().maxCoeff();
		//trama=[trama zeros(1,Lp2-L)].*w;
		// ra = real(ifft(abs(fft(trama))  . ^ 2));  --> combine these two statements in MATLAB together
		auto ftemp = fft(concat(trama, Eigen::TRowVectorX::Zero(Lp2 - L)).cwiseProduct(w).eval()).array().abs2().matrix().cast<std::complex<Eigen::TFloat>>().eval();
		auto ra_ = ifft(ftemp).real().eval();
		auto ra = (1 / ra_(0)) * ra_.head(L2);

		//  rx=ra./rw;  (MATLAB)
		//% for each element of rx, compare it with 0.0 and take the larger one
		// rx = max(rx, 0.0);
		rx = ra.cwiseQuotient(rw).cwiseMax(0.0);
		rk.setZero();
		fk.setZero();
		
		for (int j = lagmin + 1; j <= lagmax + 1; j++)
		{
			if (rx(j - 1) > 0.5*voith && rx(j - 2) < rx(j - 1) && rx(j) < rx(j - 1))
			{
				auto tmax = 0.5*(rx(j - 2) - rx(j)) / (rx(j - 2) - 2.0 * rx(j-1) + rx(j));
				auto rmax = rx(j - 1) + 0.5 * tmax * (rx(j) - rx(j - 2) + tmax * (rx(j - 2) - 2 * rx(j - 1) + rx(j)));
				if (rmax > 1)
					rmax = 1 / rmax;
				tmax = (tmax + j - 1) / fs;
				Eigen::TFloat fmax = 1 / tmax;
				rmax = rmax - octcost * std::log2(f0min * tmax);
				Eigen::Index jj = -1;
				auto rmin = rk.minCoeff(&jj);
				if (rmax > rmin)
				{
					// Here jj is the index of the min element in rk and is based on 0 starting.
					rk(jj) = rmax;
					fk(jj) = fmax;
				}
			}
		}

		Eigen::TFloat ruv = voith + std::max<Eigen::TFloat>(Eigen::TFloat(0.0), Eigen::TFloat(2.0) - std::min(Eigen::TFloat(1.0), xlmax / xgmax) / (silth / (1.0 + voith)));
		auto jj = fk.array() > 0;
		dat[k - 1].r = concat(indexByLogical(rk, jj) + octcost * (f0min / indexByLogical(fk, jj).array()).unaryExpr([](const auto& e) {return std::log2(e);}).matrix(), ruv);
		dat[k - 1].f = concat(indexByLogical(fk, jj), 0);
		dat[k - 1].ac = Eigen::TRowVectorX::Zero(dat[k - 1].f.size());
		dat[k - 1].ant = Eigen::RowVectorXi::Zero(dat[k - 1].f.size());
	}

	dat[0].ac = -dat[0].r;
	for (int k = 2; k <= pms.size(); k++)
	{
		for (int j = 1; j <= dat[k - 1].f.size(); j++)
		{
			Eigen::TFloat mincost = std::numeric_limits<Eigen::TFloat>::infinity();
			int jjmincost = -1;
			for (int jj = 1; jj <= dat[k - 2].f.size(); jj++)
			{
				Eigen::TFloat cost = dat[k - 2].ac(jj - 1) - dat[k - 1].r(j - 1);
				if (almostEqual(dat[k - 1].f(j - 1), 0) && almostEqual(dat[k - 2].f(jj - 1), 0))
					cost = cost + uvuvcost;
				else if (almostEqual(dat[k - 1].f(j - 1), 0))
					cost = cost + vuvcost;
				else if (almostEqual(dat[k - 2].f(jj - 1), 0))
					cost = cost + uvvcost;
				else
					cost = cost + octjump * std::abs(std::log2(dat[k - 1].f(j - 1) / dat[k - 2].f(jj - 1)));
				if (cost < mincost)
				{
					mincost = cost;
					jjmincost = jj;
				}
			}
			dat[k - 1].ac(j - 1) = mincost;
			dat[k - 1].ant(j - 1) = jjmincost;
		}
	}

	Eigen::TRowVectorX f0s(pms.size());
	f0s.setZero();
	// [mincost,jjmincost]=min(dat(length(dat)).ac);
	Eigen::Index jjmincost = -1;
	dat.back().ac.minCoeff(&jjmincost);
	// since jjmincost is acquired through Eigen, it is already 0-based indexing.
	f0s(f0s.size() - 1) = dat.back().f(jjmincost);
	jjmincost = dat.back().ant(jjmincost); // ant contains 1-based indices.Thereafter, jjmincost is 0-based
	for (auto k = pms.size() - 1; k >= 1; k--)
	{
		f0s(k - 1) = dat[k - 1].f(jjmincost - 1);
		jjmincost = dat[k - 1].ant(jjmincost - 1);
	}

	return f0s;
}
