#include "linphaseterm.h"
#include "sumcos.h"
#include "seq.h"

Eigen::TFloat linphaseterm(Eigen::Ref<const Eigen::TRowVectorX> aa, Eigen::Ref<const Eigen::TRowVectorX> pp, int fmaxopt)
{
	Eigen::TFloat fmax = 5000.0;
	int K = (int)std::ceil(aa.size() * fmaxopt / fmax);
	// aa_=aa(1:K); pp_=pp(1:K); MATLAB
	// Since in the remaining, aa_ and pp_ are only read and no need to create a temporary variable 
	// We just reference the subarea comprising the first K elements
	const auto aa_ = aa.head(K);
	const auto pp_ = pp.head(K);
	// Here, aa_ is actually referencing the data of aa instead of copying
	// This is more efficient. Use const to show aa_ can only be read.
	Eigen::TRowVectorX alfa = sumcos(seq(1, K).cwiseProduct(aa_), pp_, Modo::sin);

	if (alfa.size() > 1)
	{
		// combine multiple MATLAB statements here
		// In MATLAB, first compute s(k) from alfa(k), this is actually like a transformation of alpha coeffcient-wise
		// second, find the largest element in s as smax and its index imax
		int imax = 0;
		alfa.unaryExpr([&aa_, &pp_](Eigen::TFloat element) {
			return aa_.dot((seq(1, (double)pp_.size()) * element + pp_).array().cos().matrix());
		}).maxCoeff(&imax);
		return alfa(imax);
	}
	else if (alfa.size() == 1)
		return alfa.value();
	else
	{
		return 0;
		//throw std::invalid_argument("No hay alfa!!");
	}
	
	
}

