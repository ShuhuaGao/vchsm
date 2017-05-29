#include "lsfadap.h"
#include "seq.h"
#include "poly.h"
#include "concat.h"


Eigen::TRowVectorX lsfadap(Eigen::Ref<const Eigen::TVectorX> lsf)
{
	int p = (int)lsf.size();
	const auto i1 = std::complex<double>(0, 1);
	auto s1 = seq<Eigen::RowVectorXi>(2, 2, p);
	auto s2 = seq<Eigen::RowVectorXi>(1, 2, p);
	Eigen::TRowVectorXc part1(s1.size()), part2(s2.size()); // preallocate memory 
	for (Eigen::Index j = 0; j < s1.size(); j++)
	{
		part1(j) = std::exp(i1 * lsf(s1(j) - 1));
	}
	for (Eigen::Index j = 0; j < s2.size(); j++)
	{
		part2(j) = std::exp(i1 * lsf(s2(j) - 1));
	}
	Eigen::TRowVectorXc temp1(part1.size() * 2 + 1); // [1 part1 conj(part1)]
	Eigen::TRowVectorXc temp2(part2.size() * 2 + 1); // [part2 -1 conj(part2)]
	temp1(0) = 1;
	temp1.segment(1, part1.size()) = part1;
	temp1.tail(part1.size()) = part1.conjugate();
	temp2.head(part2.size()) = part2;
	temp2(part2.size()) = -1;
	temp2.tail(part2.size()) = part2.conjugate();

	auto ai = (0.5 * (poly(temp1) + poly(temp2))).real().eval();
	return ai.head(ai.size() - 1);
	
}
