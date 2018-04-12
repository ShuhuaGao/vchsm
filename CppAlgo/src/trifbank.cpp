#include "trifbank.h"
#include "seq.h"

///  hz2mel = @( hz )( 1127*log(1+hz/700) );  in mfcc.m
template<typename Derived>
auto hz2mel(const Eigen::MatrixBase<Derived>& hz)
{
	return (1 + (hz / 700).array()).log() * 1127;
}

Eigen::TFloat hz2mel(Eigen::TFloat hz)
{
	return 1127 * std::log(1 + hz / 700);
}

/// mel2hz = @( mel )( 700*exp(mel/1127)-700 ) in mfcc.m
template<typename Derived>
auto mel2hz(const Eigen::MatrixBase<Derived>& mel)
{
	return 700 * (mel.array() / 1127.0).exp() - 700;
}

Eigen::TMatrixX trifbank(int M, int K, const Eigen::TRowVector2&  R, Eigen::TFloat fs)
{
	Eigen::TFloat f_min = 0;
	Eigen::TFloat f_low = R(0);
	Eigen::TFloat f_high = R(1);
	Eigen::TFloat f_max = 0.5 * fs;
	auto f = Eigen::TVectorX::LinSpaced(K, f_min, f_max).eval();
	// h2w -> hz2mel and w2h -> mel2hz
	Eigen::TRowVectorX c = mel2hz((seq(0, M + 1).array() * ((hz2mel(f_high) - hz2mel(f_low)) / (M + 1)) + hz2mel(f_low)).matrix());
	Eigen::TMatrixX H(M, K);
	H.setZero();
	

	for (int m = 1; m <= M; m++)
	{
		for (int k = 0; k < f.size(); k++)
		{
			if (f(k) >= c(m - 1) && f(k) <= c(m)) // k -> 0-based
				H(m - 1, k) = (f(k) - c(m - 1)) / (c(m) - c(m - 1));
			else if (f(k) >= c(m) && f(k) <= c(m + 1))
				H(m - 1, k) = (c(m + 1) - f(k)) / (c(m + 1) - c(m));
		}
		
	}

	return H;
}
