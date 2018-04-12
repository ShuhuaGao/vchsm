#include "DynamicTimeWarping.h"
#include "pdist2squared.h"
#include <limits>
#include <algorithm>

/**
* Find the minimum of three numbers a, b and c, which may be NaN
*@note NaN will be neglected from comparison.
*/
std::pair<Eigen::TFloat, int> min(Eigen::TFloat a, Eigen::TFloat b, Eigen::TFloat c)
{
	int minNumIndex = 0;
	auto minNum = std::numeric_limits<Eigen::TFloat>::max();
	if (!std::isnan(a) && a < minNum)
	{
		minNum = a;
		minNumIndex = 0;
	}
	if (!std::isnan(b) && b < minNum)
	{
		minNum = b;
		minNumIndex = 1;
	}
	if (!std::isnan(c) && c < minNum)
	{
		minNum = c;
		minNumIndex = 2;
	}
	return std::make_pair(minNum, minNumIndex);
}

std::pair<std::vector<int>, std::vector<int>> DynamicTimeWarping(Eigen::Ref<const Eigen::TMatrixX> A, Eigen::Ref<const Eigen::TMatrixX> B)
{
	auto M = pdist2Squared(A, B);
	auto r = M.rows();
	auto c = M.cols();

	auto NaN = std::numeric_limits<Eigen::TFloat>::quiet_NaN();
	Eigen::TMatrixX D = Eigen::TMatrixX::Zero(r + 1, c + 1);
	D.row(0).setConstant(NaN);
	D.col(0).setConstant(NaN);
	D(0, 0) = 0;
	D.block(1, 1, r, c) = M;
	Eigen::MatrixXi phi(r, c);
	phi.setZero();

	Eigen::TFloat dmax;
	int tb;
	for (int i = 1; i <= r; i++)
		for (int j = 1; j <= c; j++)
		{
			std::tie(dmax, tb) = min(D(i - 1, j - 1), D(i - 1, j), D(i, j - 1));
			D(i, j) += dmax;
			phi(i - 1, j - 1) = tb + 1; // note: C++/Eigen is 0-based (tb)
		}

	auto i = int(r);
	auto j = int(c);
	std::vector<int> p, q;
	p.push_back(i);
	q.push_back(j);
	while (i > 1 && j > 1)
	{
		auto tb = phi(i - 1, j - 1);
		if (tb == 1)
		{
			i--;
			j--;
		}
		else if (tb == 2)
		{
			i--;
		}
		else if (tb == 3)
		{
			j--;
		}
		p.push_back(i);
		q.push_back(j);
	}
	std::reverse(p.begin(), p.end());
	std::reverse(q.begin(), q.end());
	return std::make_pair(std::move(p), std::move(q));
}
