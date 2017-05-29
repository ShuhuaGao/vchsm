#include "comp_delta.h"
#include "seq.h"

Eigen::TMatrixX comp_delta(const Eigen::Ref<const Eigen::TMatrixX>& static_coef, int DELTAWINDOW)
{
	auto N_vec = static_coef.rows();
	auto M_vec = static_coef.cols();
	// sc is the processed static_coef in matlab code
	Eigen::TMatrixX sc(N_vec + 2 * DELTAWINDOW, M_vec);
	sc.topRows(DELTAWINDOW) = static_coef.row(0).replicate(DELTAWINDOW, 1);
	sc.bottomRows(DELTAWINDOW) = static_coef.row(N_vec - 1).replicate(DELTAWINDOW, 1);
	sc.block(DELTAWINDOW, 0, N_vec, M_vec) = static_coef;
	//
	Eigen::TMatrixX delta_coef(N_vec, M_vec);
	delta_coef.setZero();
	auto i = DELTAWINDOW + 1;
	for (int j = 1; j <= DELTAWINDOW; j++)
	{
		delta_coef += j * (sc.block(i + j - 1, 0, N_vec, M_vec) - sc.block(i - j - 1, 0, N_vec, M_vec));
	}

	return delta_coef / (seq(1, DELTAWINDOW).squaredNorm() * 2);
}
