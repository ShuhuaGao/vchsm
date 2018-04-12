#include "HSMptraining.h"
#include "detalsf.h"
#include "stochalsf.h"
#include "concat.h"
#include "clustkmeans.h"
#include "StructDefinitions.h"
#include "constants.h"
#include "pinv.h"
#include "automaticth2wfw.h"
#include <algorithm>

HSMModel HSMptraining(const std::pair<PicosStructArray, PicosStructArray>& corpus, int m)
{
	constexpr int p = 14;
	// In MATLAB, load(corpus);  % two variables p1v and p2v are loaded 
	const auto& p1v = corpus.first;
	const auto& p2v = corpus.second;

	auto N = p1v.size();
	Eigen::TRowVectorX lf01(p1v.size());
	std::transform(p1v.cbegin(), p1v.cend(), lf01.data(), [](const PicosElement& e) {return std::log10(e.f0);});
	auto lf01u = (1.0 / p1v.size()) * lf01.sum();
	auto lf01v = std::sqrt(1.0 / p1v.size() * (lf01.array() - lf01u).matrix().squaredNorm());
	Eigen::TRowVectorX lf02(p2v.size());
	std::transform(p2v.cbegin(), p2v.cend(), lf02.data(), [](const PicosElement& e) {return std::log10(e.f0);});
	auto lf02u = (1.0 / p2v.size()) * lf02.sum();
	auto lf02v = std::sqrt(1.0 / p2v.size() * (lf02.array() - lf02u).matrix().squaredNorm());
	Eigen::TRowVector4 f0f12;
	f0f12 << lf01u, lf01v, lf02u, lf02v;

	Eigen::TMatrixX X, Y, Xe, Ye;
		X = detalsf(p1v, p);
		Y = detalsf(p2v, p);
		Xe = stochalsf(p1v);
		Ye = stochalsf(p2v);
	
	

	auto Z = concat<Eigen::TMatrixX>(X, Y, Eigen::Vertical); // Z =[X; Y]
	// sum along each row vector to get a column vector
	auto ugral = ((1.0 / N) * Z.rowwise().sum()).eval(); // since ugral will be used multiple times in the next, we first evaluate it.
	// some simplifications are made in the following, please check the MATLAB comments.
	auto Egral = (1.0 / N) * (Z.colwise() - ugral) * (Z.colwise() - ugral).transpose();
	auto umb = 0.02 * Egral.diagonal();
	auto pert = umb.minCoeff();

	Eigen::RowVectorXi c, c0;
	std::tie(c, c0) = clustkmeans(Z, m);

	ThzStructArray thz(m);
	std::vector<int> indi;
	for (int i = 1; i <= m;i++)
	{
		// find subZ = Z(:, c==i)
		indi.clear();
		for (int j = 0; j < c.size(); j++)
		{
			if (c[j] == i) // zero-based index
				indi.push_back(j);
		}
		int subN = (int)indi.size();
		Eigen::TMatrixX subZ(Z.rows(), subN);
		for (int j = 0; j < subN; j++)
			subZ.col(j) = Z.col(indi[j]); // zero-based index

		thz[i - 1].a = Eigen::TFloat(subN) / N;
		auto ui = Z.col(c0(i - 1) - 1);
		thz[i - 1].u = ui;

		auto Ei = (subZ.colwise() - ui) * (subZ.colwise() - ui).transpose();
		thz[i - 1].E = (1.0 / subN) * Ei + pert * Eigen::TMatrixX::Identity(2 * p, 2 * p);
	}

	Eigen::TMatrixX P(m, N);
	int seguir = 1;
	Eigen::TFloat Lant = 1;
	Eigen::TFloat cte = std::pow(2 * pi, -p);
	Eigen::TRowVectorX pxt(N);
	pxt.setZero();

	for (int i = 1; i <= m; i++)
	{
		// use reference to avoid copying
		const auto& ai = thz[i - 1].a;
		const auto& ui = thz[i - 1].u;
		const auto& Ei = thz[i - 1].E;
		auto dEi = 1 / std::sqrt(Ei.determinant());
		auto iEi = Ei.inverse().eval();
		for (int t = 1; t <= N; t++)
		{
			// Z(:,t)-ui
			auto temp = Z.col(t - 1) - ui;
			P(i - 1, t - 1) = ai * dEi * (-0.5 * temp.transpose() * iEi * temp).array().exp().value(); // note: temp.transpose() * iEi * temp is in fact a scalar
			pxt(t - 1) += cte * P(i - 1, t - 1);
		}
	}

	Lant = pxt.array().log().sum();
	auto Psum = P.colwise().sum(); // a row vector
	// each row vector of P will be divided by Psum
	P.array() = P.array().rowwise() / Psum.array();

	// Preallocate memory for size-known matrices
	Eigen::TVectorX ui(2 * p);
	Eigen::TMatrixX Ei(2 * p, 2 * p);
	while (seguir == 1)
	{
		for (int i = 1; i <= m; ++i)
			thz[i - 1].a = P.row(i - 1).sum() / N;
		for (int i = 1; i <= m; ++i)
		{
			ui.setZero();
			for (int t = 1; t <= N; t++)
				ui += P(i-1, t-1) * Z.col(t - 1);
			thz[i - 1].u = ui * (1 / (N * thz[i - 1].a));
		}

		for (int i = 1; i <= m; i++)
		{
			Ei.setZero();
			for (int t = 1; t <= N; t++)
				Ei += P(i - 1, t - 1) * (Z.col(t - 1) - thz[i - 1].u) * (Z.col(t - 1) - thz[i - 1].u).transpose();
			thz[i - 1].E = Ei * (1 / (N * thz[i - 1].a)) + pert * Eigen::TMatrixX::Identity(2 * p, 2 * p);
		}

		pxt.setZero();
		for (int i = 1; i <= m; i++)
		{
			const auto& ai = thz[i - 1].a;
			const auto& ui = thz[i - 1].u;
			const auto& Ei = thz[i - 1].E;
			auto dEi = 1 / std::sqrt(Ei.determinant());
			auto iEi = Ei.inverse().eval();
			for (int t = 1; t <= N; ++t)
			{
				P(i - 1, t - 1) = ai * dEi * std::exp((-0.5 * (Z.col(t - 1) - ui).transpose() * iEi * (Z.col(t - 1) - ui)).value());
				pxt(t - 1) += cte * P(i - 1, t - 1);
			}
		}
		auto L = pxt.array().log().sum();
		auto Psum = P.colwise().sum().eval();
		P.array() = P.array().rowwise() / Psum.array();
		if ((L - Lant) / Lant < 0.000001) // umb=0.000001; 
			seguir = 0;
		else
			Lant = L;
	}

	ThxyStructArray thxy(m);
	auto thye = thxy;
	for (int k = 1; k <= m; ++k)
	{
		thxy[k - 1].a = thz[k - 1].a;
		thxy[k - 1].u = thz[k - 1].u.head(p);
		thxy[k - 1].v = thz[k - 1].u.segment(p, p);
		thxy[k - 1].E = thz[k - 1].E.topLeftCorner(p, p);
		thxy[k - 1].R = thz[k - 1].E.block(p, 0, p, p);

		thye[k - 1].a = thz[k - 1].a;
		thye[k - 1].u = thz[k - 1].u.segment(p, p);
		thye[k - 1].v = thz[k - 1].u.head(p);
		thye[k - 1].E = thz[k - 1].E.block(p, p, p, p);
		thye[k - 1].R = thz[k - 1].E.block(0, p, p, p);
	}
	
	P = Eigen::TMatrixX::Zero(N, m);
	Eigen::TMatrixX D(N, m*p);
	D.setZero();

	for (int i = 1; i <= m; i++)
	{
		const auto& ai = thye[i - 1].a;
		const auto& ui = thye[i - 1].u;
		const auto& Ei = thye[i - 1].E;
		auto dEi = 1 / std::sqrt(Ei.determinant());
		auto iEi = Ei.inverse().eval();
		auto iEiT = iEi.transpose();
		for (int t = 1; t <= N; t++)
		{
			auto temp = (Y.col(t - 1) - ui).eval(); // column vector
			P(t - 1, i - 1) = ai * dEi * std::exp((-0.5 * temp.transpose() * iEi * temp).value());
			D.row(t - 1).segment((i - 1)*p, p) = P(t - 1, i - 1) * temp.transpose() * iEiT;
		}
	}

	auto iPsum_ = P.rowwise().sum().array().eval(); // since iPsum has already been defined in the above
	P.array().colwise() /= iPsum_.array();
	D.array().colwise() /= iPsum_.array();
	Z = (pinv(concat<Eigen::TMatrixX>(P, D, Eigen::Horizontal)) * Ye.transpose()).transpose();
	for (int i = 1; i <= m; i++)
	{
		thye[i - 1].v = Z.col(i - 1);
		thye[i - 1].R = Z.block(0, m + (i - 1)*p, Z.rows(), p);
	}

	FwxyStructArray fwxy = automaticth2wfw(thxy);

	HSMModel model;
	model.th = thxy;
	model.th2 = thye;
	model.fw = fwxy;
	model.f0f = f0f12;

	return model;
}
