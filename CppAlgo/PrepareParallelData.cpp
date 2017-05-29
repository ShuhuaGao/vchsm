#include "PrepareParallelData.h"
#include "types.h"


std::pair<PicosStructArray, PicosStructArray> PrepareParallelData(const std::vector<PicosStructArray>& sourceHSMFeatureList, const std::vector<PicosStructArray>& targetHSMFeatureList, int index)
{
	int d = 39;
	PicosStructArray p1v, p2v;


	for (std::size_t tt  = 0; tt < index; tt++)
	{
		const auto& picos_source = sourceHSMFeatureList[tt];
		auto n1 = picos_source.size();
		Eigen::TMatrixX A(n1, d);
		A.setZero();
		for (int i = 1; i <= n1; i++)
		{
			if (picos_source[i - 1].mfcc.size() != 0)
				A.row(i - 1) = picos_source[i - 1].mfcc;
		}

		const auto& picos_target = targetHSMFeatureList[tt];
		auto n2 = picos_target.size();
		Eigen::TMatrixX B(n2, d);
		B.setZero();
		for (int j = 1; j <= n2; j++)
		{
			if (picos_target[j - 1].mfcc.size() != 0)
				B.row(j - 1) = picos_target[j - 1].mfcc;
		}

		std::vector<int> p, q;
		std::tie(p, q) = DynamicTimeWarping(A, B);
		
		std::vector<int> isvv(p.size());
		for (int k = 1; k <= isvv.size(); k++)
		{
			isvv[k - 1] = (int)picos_source[p[k - 1] - 1].a.size() * (int)picos_target[q[k - 1] - 1].a.size();
		}

		//  p1v=[p1v picos_source(p(isvv))]; 
		for (int k = 1; k <= isvv.size(); k++)
		{
			if (isvv[k - 1] > 0) 
			{
				{
					p1v.push_back(picos_source[p[k - 1] - 1]);
					p2v.push_back(picos_target[q[k - 1] - 1]);
				}		
			}
		}
		
	}



	auto size = p1v.size();
	Eigen::TFloat dist = 0;
	// find elements: if dist(i) < 80, then p1v(i) and p2v(i) is chosen
	PicosStructArray p1vr, p2vr; // values to be returned
	for (std::size_t i = 1; i <= size; i++)
	{
		dist = (p1v[i - 1].mfcc - p2v[i - 1].mfcc).rowwise().norm().value();
		if (dist < 80)
		{
			p1vr.push_back(p1v[i - 1]);
			p2vr.push_back(p2v[i - 1]);
		}
	}

	return std::make_pair(std::move(p1vr), std::move(p2vr)); // save(corpus,'p1v','p2v');
}
