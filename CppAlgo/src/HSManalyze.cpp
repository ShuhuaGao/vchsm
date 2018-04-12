#include "HSManalyze.h"
#include "seq.h"
#include "f0analysisbyboersma.h"
#include "harmonicanalysis.h"
#include "stochasticanalysis.h"
#include "polarityanalysis.h"
#include "decomposephase.h"
#include "extractmfcc.h"

PicosStructArray HSManalyze(Eigen::Ref<Eigen::TRowVectorX> x, Eigen::TFloat fs)
{
	assert(fs == 16000);
	int N = 128;
	Eigen::TFloat f0min = 60;
	Eigen::TFloat f0max = 500;
	Eigen::TFloat fmax = 5000; // since fs is fixed to be 16000
	int ordenLPC = 14;
	// normalize
	x *= 0.5 / x.maxCoeff();
	int L = (int)x.size();

	auto temp = (int)std::floor((L - 3.0*fs / f0min) / N);
	Eigen::RowVectorXi pms = 1 + (int)std::ceil(1.5*fs / f0min) + N * seq<Eigen::RowVectorXi>(0, temp).array();

	auto f0s = f0analysisbyboersma(x, fs, pms, f0min, f0max);
	auto picos = harmonicanalysis(x, fs, pms, f0s, fmax);
	stochasticanalysis(x, fs, N, picos, ordenLPC);
	polarityanalysis(picos);
	decomposephase(picos);
	
	return  picos;
}
