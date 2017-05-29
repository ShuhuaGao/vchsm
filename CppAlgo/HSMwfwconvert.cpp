#include "HSMwfwconvert.h"
#include "constants.h"
#include "aaalsf.h"
#include "lsfadap.h"
#include "seq.h"
#include "concat.h"
#include "angle.h"
#include "sliceByIndices.h"
#include <algorithm>

void HSMwfwconvert(const HSMModel & model, PicosStructArray & picos)
{
	Eigen::TFloat fmax = 5000;
	const std::complex<Eigen::TFloat> i1(0, 1);

	const auto& f0f = model.f0f;
	constexpr int fs = 16000;
	auto u1 = f0f(0);
	auto v1 = f0f(1);
	auto u2 = f0f(2);
	auto v2 = f0f(3);
	Eigen::TRowVectorX f01s(picos.size());
	f01s.setZero();
	auto f02s = f01s;
	Eigen::TFloat dph = 0.0;

	for (size_t k = 1; k <= picos.size(); k++)
	{
		f01s(k - 1) = picos[k - 1].f0;
		if (f01s(k-1) > 0)
		{
			f02s(k - 1) = std::pow(10, u2 + v2 * (std::log10(f01s(k - 1)) - u1) / v1);
			if (k > 1 && f01s(k - 2) > 0)
			{
				dph = dph + (picos[k-1].pm - picos[k - 2].pm) * pi * ((f02s(k-1) - f01s(k-1)) + (f02s(k - 2) - f01s(k - 2))) / fs;
				picos[k - 1].alfa += dph;
			}
			else
			{
				dph = 0;
			}
		}
		else
		{
			f02s(k - 1) = 0;
			dph = 0;
		}
	}

	const auto& th = model.th;
	const auto& th2 = model.th2;
	auto p = th[0].u.size();
	auto m = th.size();
	auto mm = th2.size();
	
	//In C++, we cannot add a field to a struct dynamically. Just define them as variables here.
	std::vector<Eigen::TFloat> thd(th.size());
	std::vector<Eigen::TMatrixX> thI(th.size());
	std::transform(th.cbegin(), th.cend(), thd.begin(), [](const ThxyElementType& e) {return e.E.determinant();});
	std::transform(th.cbegin(), th.cend(), thI.begin(), [](const ThxyElementType& e) {return e.E.inverse();});
	std::vector<Eigen::TFloat> th2d(th2.size());
	std::vector<Eigen::TMatrixX> th2I(th2.size());
	std::transform(th2.cbegin(), th2.cend(), th2d.begin(), [](const ThxyElementType& e) {return e.E.determinant();});
	std::transform(th2.cbegin(), th2.cend(), th2I.begin(), [](const ThxyElementType& e) {return e.E.inverse();});

	const auto& fw = model.fw;
	const auto& fx = fw[0].x;
	Eigen::TRowVectorX P(m);
	Eigen::TRowVectorX PP(mm);
	auto P2 = P; // (1/Psum)*P
	for (int k = 1; k <= (int)picos.size(); k++)
	{
		if (picos[k - 1].f0 == 0)
			continue;
		auto v = aaalsf(picos[k - 1].a, picos[k - 1].f0, (int)p);
		P.setZero();
		for (size_t j = 1; j <= m; j++)
		{
			P(j - 1) = th[j - 1].a / std::sqrt(thd[j - 1]) * std::exp(-0.5 * ((v - th[j-1].u).transpose() * thI[j-1]* (v - th[j - 1].u)).value());
		}

		auto Psum = P.sum();
		P /= Psum;
		Eigen::TVectorX vt(v.size());
		vt.setZero();
		for (size_t j = 1; j <= m; j++)
		{
			vt = vt + P(j - 1) * (th[j - 1].v + th[j - 1].R * thI[j - 1] * (v - th[j - 1].u));
		}

		auto aivt = lsfadap(vt);
		Eigen::TMatrixXc eevt(p + 1, (int)std::ceil(fmax / f02s(k - 1)) - 1);
		eevt.setOnes();
		eevt(1, 0) = std::exp(-i1 * pi * f02s(k - 1) / fmax);
		for (Eigen::Index jj = 2; jj <= eevt.cols(); jj++)
		{
			eevt(1, jj - 1) = eevt(1, jj - 2) * eevt(1, 0);
		}
		for (Eigen::Index jj = 3; jj <= p + 1; jj++)
		{
			eevt.row(jj - 1) = eevt.row(jj - 2).cwiseProduct(eevt.row(1));
		}

		
		auto Afvt = (1 / (aivt * eevt).array()).matrix().eval();
		auto E1 = picos[k - 1].a.squaredNorm();
		Afvt = Afvt * std::sqrt(E1 / Afvt.dot(Afvt).real());
		auto Evt = Afvt.cwiseProduct(Afvt.conjugate()).real();
		auto aavt = Evt.cwiseSqrt();

		PP.setZero();
		for (std::size_t j = 1; j <= mm; j++)
		{
			PP(j-1) = th2[j-1].a / std::sqrt(th2d[j-1]) * std::exp(-0.5 * ((vt - th2[j - 1].u).transpose() * th2I[j - 1] * (vt - th2[j - 1].u)).value());
		}
		PP = PP / PP.sum();
		Eigen::TVectorX vtt(vt.size());
		vtt.setZero();
		for (size_t j = 1; j <= mm; j++)
		{
			vtt = vtt+ PP(j - 1) * (th2[j - 1].v + th2[j - 1].R * th2I[j - 1] * (vt - th2[j - 1].u));
		}

		Eigen::TVectorX vttaux(vtt.size() + 2);
		vttaux(0) = 0;
		vttaux(vttaux.size() - 1) = pi;
		vttaux.segment(1, vtt.size()) = vtt;
		for (Eigen::Index j = 2; j <= vttaux.size() - 1; j++)
		{
			if (vttaux(j - 1) <= vttaux(j - 2) && vttaux(j - 1) < vttaux(j) && vttaux(j - 2) < vttaux(j))
				vtt(j - 2) = 0.99 * vttaux(j - 2) + 0.01 * vttaux(j);
		}

		auto aivtt = lsfadap(vtt);
		aivtt *= (picos[k - 1].e(0) / aivtt(0));

		Eigen::TRowVectorX fy(fw[0].y.size());
		fy.setZero();
		for (std::size_t j = 1; j <= m; j++)
		{
			fy += P(j - 1) * fw[j - 1].y;
		}

		Eigen::TRowVectorX ff1(picos[k - 1].a.size() + 2);
		ff1(0) = 0.0;
		ff1(ff1.size() - 1) = fmax;
		ff1.segment(1, picos[k - 1].a.size()) = seq(1, picos[k - 1].a.size()) * f01s(k - 1);
		Eigen::TRowVectorXc ap1(picos[k - 1].a.size() + 2);
		ap1.segment(1, picos[k - 1].a.size()) = picos[k - 1].a.array() * (i1 * picos[k - 1].p).array().exp();
		ap1(0) = picos[k - 1].a(0);
		ap1(ap1.size() - 1) = 0;
		Eigen::TRowVectorX aa1(picos[k - 1].a.size() + 2);
		aa1(0) = std::log(picos[k - 1].a(0));
		aa1(aa1.size() - 1) = std::log(picos[k - 1].a.minCoeff());
		aa1.segment(1, picos[k - 1].a.size()) = picos[k - 1].a.array().log();
		Eigen::TRowVectorX aa2((int)std::ceil(fmax / f02s(k - 1)) - 1);
		aa2.setZero();
		auto pp2 = aa2;
		int jjant1 = 1;
		int jjant2 = 1;
		Eigen::TFloat scale = 0.0;
		for (Eigen::Index j = 1; j <= aa2.size(); j++)
		{
			auto f2j = j * f02s(k - 1);
			Eigen::Index jj = 0;
			for (jj = jjant1; jj <= fy.size() - 1; jj++)
			{
				if (f2j >= fy(jj - 1) && f2j < fy(jj))
				{
					jjant1 = (int)jj;
					break;
				}
			}
			auto f1j = fx(jj - 1) + (fx(jj) - fx(jj-1))*(f2j - fy(jj-1)) / (fy(jj) - fy(jj-1));
			if (f1j < fmax)
			{
				for (jj = jjant2; jj <= ff1.size() - 1; jj++)
				{
					if (f1j >= ff1(jj - 1) && f1j < ff1(jj))
					{
						jjant2 = (int)jj;
						break;
					}
				}
				aa2(j-1) = std::exp(aa1(jj-1) + (aa1(jj) - aa1(jj-1)) * (f1j - ff1(jj-1)) / (ff1(jj) - ff1(jj-1)));
				pp2(j-1) = angle(ap1(jj - 1) + (ap1(jj) - ap1(jj-1))*(f1j - ff1(jj-1)) / (ff1(jj) - ff1(jj-1)));
			}
			else
			{
				if (scale == 0)
				{
					scale = aa2(j - 2) / std::abs(Afvt(j - 2));
				}
				aa2(j - 1) = scale * std::abs(Afvt(j - 1));
				pp2(j - 1) = angle(Afvt(j - 1));
			}
		}

		// modo is always 'fwa' in the MATLAB program
		Eigen::TFloat flimsm = 1000.0;
		int nf = (int)std::ceil(flimsm / f02s(k - 1) - 1);
		auto seqnf = seq<Eigen::RowVectorXi>(-nf, nf).eval();
		Eigen::TRowVectorX win = 1 - (f02s(k - 1) / flimsm) * seqnf.cast<Eigen::TFloat>().array();
		win /= win.sum();
		auto g = (aavt.array() / aa2.array()).log().matrix().eval();
		auto gsm = g;
		Eigen::RowVectorXi ind(seqnf.size());
		for (Eigen::Index j = 1; j <= g.size(); j++)
		{
			ind = seqnf.array() + j;
			if (j <= nf)
			{
				for (Eigen::Index jj = 1; jj <= ind.size(); ++jj)
				{
					if (ind(jj - 1) < 0)
						ind(jj - 1) = -ind(jj - 1);
					else if (ind(jj - 1) == 0)
						ind(jj - 1) = 1;
					else
						break;
				}
			}
			else if (j > g.size() - nf)
			{
				for (auto jj = ind.size(); jj >= 1; jj--)
				{
					if (ind(jj - 1) > g.size() + 1)
						ind(jj - 1) = 2 * (int)g.size() + 2 - ind(jj - 1);
					else if (ind(jj - 1) == g.size() + 1)
						ind(jj - 1) = (int)g.size();
					else
						break;
				}
			}
			gsm(j - 1) = sliceByIndices(g, ind, IndexBase::One).dot(win);
		}
		g = gsm.array().exp();
		aa2 = aa2.array() * g.array();
		aa2 *= std::sqrt(E1 / aa2.squaredNorm());

		picos[k - 1].f0 = f02s(k - 1);
		picos[k - 1].a = aa2;
		picos[k - 1].p = pp2;
		picos[k - 1].e = aivtt;
	}

}
