#include "bylevdurb.h"

Eigen::TRowVectorX bylevdurb(const Eigen::Ref<const Eigen::TRowVectorX>& R)
{
	// Eigen::Index --> std::ptrdiff_t --> int
	Eigen::Index ordenLPC = R.cols() - 1;
	Eigen::TRowVectorX ai = Eigen::TRowVectorX::Zero(ordenLPC);
	Eigen::TFloat e = R(0);

	for (Eigen::Index u = 1; u <= ordenLPC; u++)
	{
		Eigen::TFloat k = 0;
		if (u == 1)
		{
			k = R(1) / R(0);
			ai(0) = k;
			e = (1 - k * k) * e;
		}
		else
		{
			//k = R(u+1) - dot(ai(1:1:u-1), R(u:-1:2)); (MATLAB)
			k = R(u) - ai.head(u - 1).dot(R.segment(1, u - 1).reverse());
			k = k / e;
			ai(u - 1) = k;
			// ai appears on both sides of assignment. eval() is added to avoid aliasing. 
			ai.head(u - 1) = ai.head(u - 1) - k * ai.head(u - 1).reverse().eval();		
			e = (1 - k * k) * e;
		}
	}

	return  (Eigen::TRowVectorX(ai.cols() + 1) << 1, -ai).finished() / sqrt(std::abs(R(0) - ai.dot(R.segment(1, ordenLPC))));
}
