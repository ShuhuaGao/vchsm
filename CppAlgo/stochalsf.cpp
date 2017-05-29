#include "stochalsf.h"
#include "concat.h"
#include "angle.h"
#include "roots.h"
#include <algorithm>

Eigen::TMatrixX stochalsf(const PicosStructArray& picos)
{
	auto p = picos[0].e.size() - 1;
	Eigen::TMatrixX LSF(p, picos.size());
	LSF.setZero();

	for (int k = 1; k <= picos.size(); k++)
	{
		/*auto ai = picos[k - 1].e;
		ai *= 1 / ai(0);*/
		auto az1 = concat(picos[k-1].e / picos[k-1].e(0), 0);
		auto az2 = az1.reverse();
		// lsf=angle([roots(az1+az2); roots(az1-az2)]); 
		Eigen::TVectorX lsf = angle(concat<Eigen::TVectorXc>(roots(az1 + az2), roots(az1 - az2), Eigen::Vertical));
		std::sort(lsf.data(), lsf.data() + lsf.size()); // sort lsf in ascending order
		LSF.col(k - 1) = lsf.segment(lsf.size() - p - 1, p);
	}
	return LSF;
}
