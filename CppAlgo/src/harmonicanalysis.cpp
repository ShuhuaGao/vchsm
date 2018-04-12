#include "harmonicanalysis.h"
#include "constants.h"
#include "angle.h"
#include "seq.h"


PicosStructArray harmonicanalysis(Eigen::Ref<const Eigen::TRowVectorX>  x, Eigen::TFloat fs, Eigen::Ref<const Eigen::RowVectorXi> pms, Eigen::Ref<const Eigen::TRowVectorX> f0s, Eigen::TFloat fmax)
{
	// define a struct array, size: the length of vector pms
	PicosStructArray picos(pms.size());
	auto size = pms.size();

	for (int k = 1; k <= size; k++)
	{
		picos[k - 1].pm = pms(k - 1);
		picos[k - 1].f0 = f0s(k - 1);
		if (f0s(k - 1) > 0)
		{
			int Lw = (int)std::ceil(2.2 * fs / std::min(f0s(k - 1), Eigen::TFloat(150.0)));
			int Lw2 = (int)std::floor(Lw / 2.0); 
			Lw = 2 * Lw2 + 1;
			// trama2T=x(pms(k)+(-Lw2:Lw2)); MATLAB
			// equivalent to x(pms(k)-Lw2 : pms(k) + Lw2)
			auto trama2T = x.segment(pms(k - 1) - Lw2 - 1, 2 * Lw2 + 1);
			// classical Hamming window
			// win=transpose(0.54+0.46*cos(pi*(-Lw2:Lw2)/Lw2));
			// no direct support in Eigen to generate vector gv = -Lw2:Lw2, use LinSpace instead
			const Eigen::TRowVectorX gv = seq(-Lw2, Lw2);
			Eigen::TVectorX win = ((gv * pi / Lw2).array().cos() * 0.46 + 0.54).transpose();

			int K = (int)std::ceil(fmax / f0s(k - 1)) - 1;
			// since h will be double columned in the next, we first allocate it here for efficiency
			Eigen::TMatrixXc h(Lw, 2*K); // h is a matrix of complex Eigen::TFloat
			h.setZero();
			// first column
			std::complex<Eigen::TFloat> i1(0, 1); // 0 + 1i
			// h(:,1)=transpose(   exp( 1i*2*pi*f0s(k)*(-Lw2:Lw2)/fs )  );  in MATLAB, note transpose() is purely transpose with no conjugate 
			h.col(0) = (gv / fs * i1 * 2 * pi * f0s(k - 1)).array().exp().transpose();

			for (int kk = 2; kk <= K; kk++)
			{
				h.col(kk - 1) = h.col(kk - 2).cwiseProduct(h.col(0));
			}

			for (int kk = 1; kk <= K; kk++)
			{
				h.col(kk - 1) = h.col(kk - 1).cwiseProduct(win);
			}
			// h=[h conj(h)]; MATLAB
			h.rightCols(K) = h.leftCols(K).conjugate();
			// coef=(h'*h)\(h'*trama2T);  MATLAB, note: since h is complex matrix, h' is conjugate transpose in MATLAB 
			// use adjoint to find the conjugate transpose in Eigen
			// least square solution 
			auto t1 = h.adjoint().eval();
			auto t2 = (t1 * h).eval();
			auto t3 = (t1 * trama2T.transpose().cwiseProduct(win)).eval();
			auto coef = t2.llt().solve(t3).eval(); // lu() gives a temporary, must eval()
			auto coef1K = coef.head(K);
			picos[k - 1].a = 2 * coef1K.cwiseAbs().transpose();
			picos[k - 1].p = angle(coef1K).transpose();
			
		}
	}
	return picos;
}
