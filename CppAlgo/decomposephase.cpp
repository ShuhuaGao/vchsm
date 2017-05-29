#include "decomposephase.h"
#include "seq.h"
#include "linphaseterm.h"

void decomposephase(PicosStructArray& picos)
{
	int fmaxopt = 1000;

	// Like MATLAB, iterate over picos
	for (int i = 0; i < picos.size(); i++)
	{
		auto& pk = picos[i];
		if (pk.f0 > 0)
		{
			Eigen::TFloat alfa = linphaseterm(pk.a, pk.p, fmaxopt);
			pk.alfa = -alfa;
			pk.p +=  seq(1, (int)pk.p.size()) * alfa;
		}
	}
}
