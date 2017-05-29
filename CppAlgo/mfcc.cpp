#include "mfcc.h"
#include "pow2roundup.h"
#include "seq.h"
#include "constants.h"
#include "filter.h"
#include "vec2frames.h"
#include "FFT.h"
#include "trifbank.h"

#ifdef MSVC
#pragma warning(disable : 4503)
#endif

Eigen::TMatrixX mfcc(Eigen::Ref<Eigen::TVectorX> speech, Eigen::TFloat fs)
{
	Eigen::TFloat Tw = 25;
	Eigen::TFloat Ts = 8;
	Eigen::TFloat alpha = 0.97;
	Eigen::Matrix<Eigen::TFloat, 1, 2> R;
	R << 300, 3700;
	int M = 20;
	int N = 13;
	int L = 22;

	// Explode samples to the range of 16 bit shorts
	if (speech.cwiseAbs().maxCoeff() <= 1)
		speech *= Eigen::TFloat(1 << 15);


	int Nw = static_cast<int>(std::round(1e-3 * Tw * fs));
	int Ns = static_cast<int>(std::round(1e-3 * Ts * fs));

	int nfft = pow2Roundup(Nw);
	int K = nfft / 2 + 1;

	auto dctm = [](int N, int M) {
		auto rep1 = seq<Eigen::TVectorX>(0, N - 1).replicate(1, M);
		auto rep2 = ((seq<Eigen::TRowVectorX>(1, M).array() - 0.5) * (pi / M)).replicate(N, 1);
		return (std::sqrt(2.0 / M) * (rep1.array() * rep2.array()).cos()).eval();
	};

	auto ceplifter = [](int N, int L) {
		return (1 + (0.5 * L) * (seq(0, N - 1) * (pi / L)).array().sin()).eval();
	};

	filter(alpha, speech);

	auto frames = vec2frames(speech, Nw, Ns);
	auto MAG = fft(frames, nfft, 1).cwiseAbs().eval();
	auto H = trifbank(M, K, R, fs);
	// analyze the matrix product here
	auto FBE = H * MAG.topRows(K);
	auto DCT = dctm(N, M);
	auto CC_ = DCT.matrix() * FBE.array().log().matrix();
	auto lifter = ceplifter(N, L);
	auto CC = lifter.matrix().asDiagonal() * CC_.matrix();

	return CC;
}