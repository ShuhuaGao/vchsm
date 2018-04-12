#include "extractmfcc.h"
#include "mfcc.h"
#include "comp_delta.h"

void extractmff(Eigen::Ref<Eigen::TRowVectorX> x, Eigen::TFloat fs, PicosStructArray& picos)
{
	auto Npm = picos.size();
	int N = 128;
	auto trama = x.head(picos[Npm - 1].pm - static_cast<int>(std::floor(N / 2.0)) + N - 1).transpose();

	auto MFCCs = mfcc(trama, fs);
	auto delta_mfcc = comp_delta(MFCCs.transpose(), 3);
	auto delta_delta_mfcc = comp_delta(delta_mfcc, 2);
	Eigen::TMatrixX mfcc_v(MFCCs.rows() + delta_mfcc.cols() + delta_delta_mfcc.cols(), MFCCs.cols());
	mfcc_v << MFCCs, delta_mfcc.transpose(), delta_delta_mfcc.transpose();

	for (Eigen::Index k = 1; k <= Eigen::Index(Npm); k++)
	{
		picos[k - 1].mfcc = mfcc_v.col(k - 1);
	}
	

}
