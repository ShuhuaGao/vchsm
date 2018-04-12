#include "polarityanalysis.h"
#include "constants.h"


int polarityanalysis(PicosStructArray& picos)
{
	int pol = 0;
	Eigen::TFloat polE = 0.0;

	for (int k = 1; k <= picos.size(); k++)
	{
		if (picos[k - 1].f0 > 0)
		{
			Eigen::TFloat E = picos[k - 1].a * picos[k - 1].a.transpose();
			Eigen::TFloat alfa = 2 * picos[k - 1].p(0) - picos[k - 1].p(1);
			alfa = alfa - 2 * pi* std::floor(alfa / (2 * pi));
			// dP=min(abs([alfa 2*pi-alfa]));  MATLAB
			Eigen::TFloat dP = std::min(std::abs(alfa), std::abs(2 * pi - alfa));
			Eigen::TFloat dN = std::abs(alfa - pi);
			if (dN < dP)
			{
				pol = pol - 1;
				polE = polE - E;
			}
			else if (dP < dN)
			{
				pol = pol + 1;
				polE = polE + E;
			}
			
		}
	}

	if (pol < 0)
		pol = -1;
	else if (pol > 0)
		pol = 1;
	else if (polE < 0)
		pol = -1;
	else
		pol = 1;

	if (pol == -1)
	{
		// like MATLAB, enumerate each element in picos 
		for (PicosElement& element : picos)
		{
			if (element.f0 > 0)
				element.p.array() += pi;
		}
	}
	return pol;
}