#include "clustkmeans.h"
#include "seq.h"
#include "StructDefinitions.h"
#include <limits>
#include <vector>



std::pair<Eigen::RowVectorXi, Eigen::RowVectorXi> clustkmeans(Eigen::Ref<const Eigen::TMatrixX> X, int n)
{
	auto N = X.cols();
	Eigen::RowVectorXi c0 = seq<Eigen::RowVectorXi>(1, n) * int(std::floor(N / Eigen::TFloat(n)));
	Eigen::RowVectorXi c (N);
	c.setZero();
	Eigen::RowVectorXi Nj (n);
	Nj.setZero();

	for (int k = 1; k <= N; k++)
	{
		auto dmin = std::numeric_limits<Eigen::TMatrixX::Scalar>::infinity();
		int jmin = 0;
		for (int j = 1; j <= n; j++)
		{
			auto dact = (X.col(k - 1) - X.col(c0(j - 1) - 1)).squaredNorm();
			if (dact < dmin)
			{
				dmin = dact;
				jmin = j;
			}
		}
		c(k - 1) = jmin;
		Nj(jmin - 1)++;
	}

	auto copiac = c;
	int seguir = 1;
	int iter = 0;
	int maxiters = 200;
	std::vector<AuxElementType> aux(n);
	while (seguir == 1)
	{
		for (int k = 1; k <= n; k++)
		{
			aux[k - 1].act = 1;
			aux[k - 1].ind = Eigen::RowVectorXi::Zero(Nj(k - 1));
			aux[k - 1].sumd = Eigen::TRowVectorX::Zero(Nj(k - 1));
		}

		for (int k = 1; k <= N; k++)
		{
			auto cls = c(k - 1);
			auto act = aux[cls - 1].act;
			aux[cls - 1].ind(act - 1) = k;
			for (int j = 1; j <= act - 1; j++)
			{
				auto dact = (X.col(k - 1) - X.col(aux[cls - 1].ind(j - 1))).squaredNorm();
				aux[cls - 1].sumd(j - 1) += dact;
				aux[cls - 1].sumd(act - 1) += dact;
			}
		}

		for (int k = 1; k <= n; k++)
		{
			int j = 0;
			aux[k - 1].sumd.minCoeff(&j);
			c0(k - 1) = aux[k - 1].ind(j);
		}

		Nj.setZero();
		for (int k = 1; k <= N; k++)
		{
			auto dmin = std::numeric_limits<Eigen::TMatrixX::Scalar>::infinity();
			int jmin = 0;
			for (int j = 1; j <= n; j++)
			{
				auto dact = (X.col(k - 1) - X.col(c0(j - 1) - 1)).squaredNorm();
				if (dact < dmin)
				{
					dmin = dact;
					jmin = j;
				}
			}
			c(k - 1) = jmin;
			Nj(jmin - 1)++;
		}

		if ((c.array() == copiac.array()).all())
			seguir = 0;
		else
		{
			copiac = c;
		}

		iter++;
		if (iter >= maxiters)
			seguir = 0;
	}


	return std::make_pair(c, c0);
}
