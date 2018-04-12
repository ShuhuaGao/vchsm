#include "automaticth2wfw.h"
#include "poly.h"
#include "concat.h"
#include "seq.h"
#include "roots.h"
#include "angle.h"
#include "constants.h"
#include "sliceByIndices.h"
#include "StructDefinitions.h"
#include <tuple>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>


/**
R=roots(ai).'; % nonconjugate transpose

% phase and real value for the poles
Ra=angle(R); % rad frequency
Rm=abs(R);

ind=find((Ra>0.0)&(Ra<pi));
Ra=Ra(ind);
Rm=Rm(ind);
R=R(ind);

% store the poles in increasing order
[Ra,ind]=sort(Ra);
Rm=Rm(ind);
R=R(ind);

% delet the first pole
if length(R)>=p/2
ind=2:length(R);
Ra=Ra(ind);
Rm=Rm(ind);
R=R(ind);
end;
*/
// return Ra, Rm and R
std::tuple<Eigen::TRowVectorX,  Eigen::TRowVectorXc> get_Ra__R(const Eigen::TRowVectorX& ai, int p)
{
	auto R = roots(ai).transpose().eval();
	auto Ra = angle(R).eval();
	// find Ra >0 and Ra < pi, then sort it. combine these two operations.
	// 1. find the index of Ra elements that are in (0, pi)
	std::vector<int> ind(Ra.size());
	std::iota(ind.begin(), ind.end(), 0);
	auto endIt = std::copy_if(ind.cbegin(), ind.cend(), ind.begin(), [&Ra](int index) {return Ra(index) > 0 && Ra(index) < pi;});
	// 2. sort these found elements in ascending order
	std::sort(ind.begin(), endIt, [&Ra](int ind1, int ind2) {return Ra(ind1) < Ra(ind2);});
	// 3. % delet the first pole if length(R)>=p/2 
	auto beginIt = ind.begin();
	if (std::distance(beginIt, endIt) >= p / 2)
		++beginIt;
	// 4. now the range [beginIt, endIt) contains the indices to be chosen for Ra, Rm, and R
	int size = (int)std::distance(beginIt, endIt);
	Eigen::TRowVectorX   Ra_(size);
	Eigen::TRowVectorXc R_(size);
	int i = 0;
	for (auto it = beginIt; it != endIt; ++it)
	{
		R_(i) = R(*it);
		Ra_(i) = Ra(*it);
		++i;
	}
	return std::make_tuple(Ra_,  R_);
}

// get ai, cc and R
std::tuple<Eigen::TRowVectorX, Eigen::TRowVectorX, Eigen::TRowVectorXc> get_ai_cc_R(Eigen::Ref<const Eigen::TRowVectorX> lsf, int p)
{
	constexpr std::complex<Eigen::TFloat> i1(0, 1);
	Eigen::Map<const Eigen::TRowVectorX, 0, Eigen::InnerStride<2>> lsf2p(lsf.data() + 1, 7); // (2,4,6,8,10,12,14) index in MATLAB
	Eigen::Map<const Eigen::TRowVectorX, 0, Eigen::InnerStride<2>> lsf1p(lsf.data(), 7); // (1,3,5,7,9,11,13) index in MATLAB
	auto part1 = (lsf2p * i1).array().exp().matrix().eval();
	auto part2 = (lsf1p * i1).array().exp().matrix().eval();

	auto pSize = part1.size();
	Eigen::TRowVectorXc roots(2 * pSize + 1);
	// 1st polynomial
	roots(0) = std::complex<Eigen::TFloat>(1, 0);
	roots.segment(1, pSize) = part1;
	roots.tail(pSize) = part1.conjugate();
	auto polynomial1 = poly(roots).real().eval();
	// 2nd polynomial
	roots.head(pSize) = part2;
	roots(pSize) = std::complex<Eigen::TFloat>(-1, 0);
	roots.tail(pSize) = part2.conjugate();
	auto polynomial2 = poly(roots).real().eval();
	auto ai_temp = 0.5 * (polynomial1 + polynomial2);
	auto ai =  ai_temp.head(ai_temp.size() - 1).eval();

	Eigen::TRowVectorX cc(lsf.size());
	cc.setZero();
	auto ai_ = -ai.segment(1, p);
	cc(0) = ai_(0);
	for (int n = 2; n <= p; n++)
	{
		cc(n - 1) = ai_(n - 1) + ((1 - seq(1, n - 1).array() / n) * ai_.head(n - 1).array() * cc.head(n - 1).reverse().array()).sum();
	}

	Eigen::TRowVectorX Ra;
	Eigen::TRowVectorXc R;
	std::tie(Ra,  R) = get_Ra__R(ai, p);

	Eigen::TMatrixX subais(R.size(), 3);
	subais.col(0) = Eigen::TVectorX(R.size());
	subais.col(0).setOnes();
	subais.col(1) = -2 * R.transpose().real();
	subais.col(2) = (R.array() * R.conjugate().array()).real(); // in fact, R .* R' is already a real vector
	Eigen::TMatrixXc emiwj(3, R.size());
	emiwj.setOnes();
	emiwj.row(1) = (-i1 * Ra).array().exp();
	emiwj.row(2) = emiwj.row(1).cwiseProduct(emiwj.row(1));
	auto mods = 1 / (subais * emiwj).cwiseAbs().array();
	Eigen::RowVectorXi ind(mods.cols()); // col-wise maximum element index starting from 0
	for (int i = 0; i < mods.cols(); i++)
	{
		int index = -1;
		mods.col(i).maxCoeff(&index);
		ind(i) = index;
	}
	auto isEqual = (ind.array() == seq<Eigen::RowVectorXi>(0, (int)R.size() - 1).array()).eval();
	auto count = isEqual.count(); // find the number of non-zero elements (true)
	if (count < R.size())
	{
		Eigen::TRowVectorXc R_(count);
		int current = 0;
		for (int i = 0; i < isEqual.size(); i++)
		{
			if (isEqual(i))
			{
				R_(current) = R(i);
				current++;
			}
		}
		return std::make_tuple(ai, cc, R_);
	}
	return std::make_tuple(ai, cc, R);
}


FwxyStructArray automaticth2wfw(const ThxyStructArray & th, int fs2, int fmax)
{
	const auto m = th.size();
	constexpr int p = 14; // from the MATLAB code, we see p is fixed to be 14
	assert(p == th[0].u.size());
	PolosStructArray polos(m);
	constexpr std::complex<Eigen::TFloat> i1(0, 1);

	for (int k = 1; k <= m; k++)
	{
		//auto lsf = th[k - 1].u.transpose();
		std::tie(polos[k - 1].ai1, polos[k - 1].cc1, polos[k - 1].p1) = get_ai_cc_R(th[k - 1].u.transpose(), p);
		// second part, Y target vector
		//lsf = th[k - 1].v.transpose();
		std::tie(polos[k - 1].ai2, polos[k - 1].cc2, polos[k - 1].p2) = get_ai_cc_R(th[k - 1].v.transpose(), p);
	}

	const Eigen::RowVector3i V1(1, -1, 1);
	const Eigen::RowVector3i V2(-1, 1, -1);
	Eigen::TMatrixX mat(2 * p, 2 * p);
	mat.setZero();

	for (int k = 1; k <= m; k++)
	{
		int Np1 = (int)polos[k - 1].p1.size();
		int Np2 = (int)polos[k - 1].p2.size();
		const auto& cc1 = polos[k - 1].cc1;
		const auto& cc2 = polos[k - 1].cc2;

		auto ccd1 = cc1.cwiseProduct(seq(1, (int)cc1.size())).eval();
		auto ccd2 = cc2.cwiseProduct(seq(1, (int)cc2.size())).eval();
		auto E = std::sqrt(cc1.squaredNorm() / ccd1.squaredNorm());
		ccd1 *= E;
		ccd2 *= E;

		auto p1 = angle(polos[k - 1].p1).eval();
		auto p2 = angle(polos[k - 1].p2).eval();
		// compute the following wx and xy
		auto get_wx_wy = [](const Eigen::TRowVectorX& p, const Eigen::RowVectorXi& gr) {
			Eigen::TRowVectorX r(2 + gr.size());
			r(0) = 0;
			r(r.size() - 1) = pi;
			// get p(gr), here gr works as an index array
			for (int i = 0; i < gr.size(); i++)
				r(i + 1) = p(gr(i) - 1);
			return r;
		};
		
		Eigen::RowVectorXi gr1opt;
		Eigen::RowVectorXi gr2opt;
		Eigen::TFloat dmin = std::numeric_limits<Eigen::TFloat>::infinity();

		auto jBegin = std::min<int>({ 3, Np1, Np2 });
		auto jEnd = std::min(Np1, Np2);
		for (int j = jBegin; j <= jEnd; j++)
		{
			auto gr1 = seq<Eigen::RowVectorXi>(1, j).eval(); // gr1 -> 1:j
			auto gr1max = ((Np1 + 1) - gr1.reverse().array()).matrix().eval();
			while (true)
			{
				auto gr2 = seq<Eigen::RowVectorXi>(1, j).eval(); // gr2 -> 1:j
				auto gr2max = ((Np2 + 1) - gr2.reverse().array()).matrix().eval();
				while (true)
				{
					auto wx = get_wx_wy(p1, gr1);
					auto wy = get_wx_wy(p2, gr2);
					int demasiadoirregular = 0;
					if (wx.size() >= 3)
					{
						auto signos = (wx - wy).cwiseSign().cast<int>().eval();
						signos(signos.size() - 1) = -signos(signos.size() - 2);
						

						for (int jj = 1; jj <= signos.size() - 2; ++jj)
						{
							auto s = signos.segment(jj - 1, 3);
							if (s == V1 || s == V2)
							{
								demasiadoirregular = 1;
								break;
							}
						}
					}

					if (demasiadoirregular == 0)
					{
						Eigen::TFloat dact1 = 0, dact2 = 0;
						int bound = (int)wx.size() - 1;
						for (int jj = 1; jj <= bound; jj++)
						{
							auto A = (wy(jj) - wy(jj - 1)) / (wx(jj) - wx(jj - 1));
							auto B = wy(jj - 1) - A * wx(jj - 1);
							for (int j1 = 1; j1 <= 2 * p; j1++)
							{
								for (int j2 = 1; j2 <= j1; j2++)
								{
									if (j1 <= p)
									{
										if (j1 == j2)
											mat(j1-1,j2-1) = 0.5*cc1(j1-1)*cc1(j2-1)*(((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj)) + wx(jj)) - ((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj-1)) + wx(jj-1)));
										else
											mat(j1-1, j2-1) = 0.5*cc1(j1-1)*cc1(j2-1)*(((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj)) + (1.0 / (j1 - j2))*std::sin((j1 - j2)*wx(jj))) - ((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj-1)) + (1.0 / (j1 - j2))*std::sin((j1 - j2)*wx(jj-1))));
									}
									else if (j2 > p)
									{
										auto j1_ = j1 - p;
										auto j2_ = j2 - p;
										if (j1_ == j2_)
											mat(j1-1, j2-1) = 0.5*cc2(j1_-1)*cc2(j2_-1)*(((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj) + B)) + wx(jj)) - ((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj-1) + B)) + wx(jj-1)));
										else
											mat(j1-1, j2-1) = 0.5*cc2(j1_-1)*cc2(j2_-1)*(((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj) + B)) + (1.0 / ((j1_ - j2_)*A))*std::sin((j1_ - j2_)*(A*wx(jj) + B))) - ((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj-1) + B)) + (1.0 / ((j1_ - j2_)*A))*std::sin((j1_ - j2_)*(A*wx(jj-1) + B))));
									}
									else
									{
										auto j1_ = j1 - p;
										mat(j1-1, j2-1) = -0.5*cc1(j1_-1)*cc2(j2-1)*(((1.0 / (j1_ + A*j2))*std::sin((j1_ + A*j2)*wx(jj) + B*j2) + (1 / (j1_ - A*j2))*sin((j1_ - A*j2)*wx(jj) - B*j2)) - ((1.0 / (j1_ + A*j2))*std::sin((j1_ + A*j2)*wx(jj-1) + B*j2) + (1.0 / (j1_ - A*j2))*std::sin((j1_ - A*j2)*wx(jj-1) - B*j2)));
									}

									if (j1 != j2)
										mat(j2 - 1, j1 - 1) = mat(j1 - 1, j2 - 1);
								}
							}
							dact1 = dact1 + (1 + A) * mat.sum();

							for (int j1 = 1; j1 <= 2 * p; j1++)
							{
								for (int j2 = 1; j2 <= j1; j2++)
								{
									if (j1 <= p)
									{
										if (j1 == j2)
											mat(j1-1, j2-1) = -0.5*ccd1(j1-1)*ccd1(j2-1)*(((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj)) - wx(jj)) - ((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj-1)) - wx(jj-1)));
										else 
											mat(j1-1, j2-1) = -0.5*ccd1(j1-1)*ccd1(j2-1)*(((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj)) - (1.0 / (j1 - j2))*std::sin((j1 - j2)*wx(jj))) - ((1.0 / (j1 + j2))*std::sin((j1 + j2)*wx(jj-1)) - (1.0 / (j1 - j2))*std::sin((j1 - j2)*wx(jj-1))));
									}
									else if (j2 > p)
									{
										auto j1_ = j1 - p;
										auto j2_ = j2 - p;
										if (j1_ == j2_)
											mat(j1-1, j2-1) = -0.5*ccd2(j1_-1)*ccd2(j2_-1)*(((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj) + B)) - wx(jj)) - ((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj-1) + B)) - wx(jj-1)));
										else 
											mat(j1-1, j2-1) = -0.5*ccd2(j1_-1)*ccd2(j2_-1)*(((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj) + B)) - (1.0 / ((j1_ - j2_)*A))*std::sin((j1_ - j2_)*(A*wx(jj) + B))) - ((1.0 / ((j1_ + j2_)*A))*std::sin((j1_ + j2_)*(A*wx(jj-1) + B)) - (1.0 / ((j1_ - j2_)*A))*std::sin((j1_ - j2_)*(A*wx(jj-1) + B))));

									}
									else
									{
										auto j1_ = j1 - p;
										mat(j1-1, j2-1) = 0.5*ccd1(j1_-1)*ccd2(j2-1)*(((1.0 / (j1_ + A*j2))*std::sin((j1_ + A*j2)*wx(jj) + B*j2) - (1.0 / (j1_ - A*j2))*std::sin((j1_ - A*j2)*wx(jj) - B*j2)) - ((1 / (j1_ + A*j2))*std::sin((j1_ + A*j2)*wx(jj-1) + B*j2) - (1.0 / (j1_ - A*j2))*std::sin((j1_ - A*j2)*wx(jj-1) - B*j2)));
									}
									if (j1 != j2)
										mat(j2 - 1, j1 - 1) = mat(j1 - 1, j2 - 1);
								}
							}
							dact2 = dact2 + (1 + A) * mat.sum();
						}

						auto dact = dact1 + dact2;
						if (dact < dmin)
						{
							dmin = dact;
							gr1opt = gr1;
							gr2opt = gr2;
						}
					}

					if (gr2 == gr2max) // matrix equality check --> bool (true if all coefficients are equal)
						break;
					else
					{
						for (int jj = j; jj >= 1; jj--)
						{
							if (gr2(jj - 1) < gr2max(jj - 1))
							{
								gr2(jj - 1)++;
								for (int jjj = jj + 1; jjj <= j; jjj++)
								{
									gr2(jjj - 1) = gr2(jjj - 2) + 1;
								}
								break;
							}
						}
					}

				}

				if (gr1 == gr1max)
					break;
				else
				{
					for (int jj = j; jj >= 1; jj--)
					{
						if (gr1(jj - 1) < gr1max(jj - 1))
						{
							gr1(jj - 1)++;
							for (int jjj = jj + 1; jjj <= j; jjj++ )
							{
								gr1(jjj - 1) = gr1(jjj - 2) + 1;
							}
							break;
						}
					}
				}
			}

		} // end of for (int j = jBegin; j <= jEnd; j++)
		// Here we need to reserve only partial elements of p1 and p2 by the indices gr1opt and gr2opt respectively.
		polos[k - 1].p1 = sliceByIndices(polos[k - 1].p1, gr1opt, IndexBase::One);
		polos[k - 1].p2 = sliceByIndices(polos[k - 1].p2, gr2opt, IndexBase::One);
	}

	auto ff1 = seq<Eigen::RowVectorXi>(0, 100, fs2).eval();
	Eigen::TRowVectorX ff2(ff1.size());
	ff2.setZero();
	// fwxy(1:m)=struct('x',[],'y',[]); 
	FwxyStructArray fwxy(m);

	for (int k = 1; k <= m; k++)
	{
		// note: ff1r(length(ff1r)) is the last element of ff1r
		auto angle1 = (fmax / pi) * angle(polos[k - 1].p1);
		auto angle2 = (fmax / pi) * angle(polos[k - 1].p2);
		auto lastElement1 = angle1(angle1.size() - 1);
		auto lastElement2 = angle2(angle2.size() - 1);
		Eigen::TRowVectorX ff1r(angle1.size() + 2);
		Eigen::TRowVectorX ff2r(angle2.size() + 2);
		if ( lastElement1 > lastElement2 )
		{
			ff1r << 0, angle1, fs2 * lastElement1 / lastElement2;
			ff2r << 0, angle2, fs2;
		}
		else
		{
			ff2r << 0, angle2, fs2 * lastElement2 / lastElement1;
			ff1r << 0, angle1, fs2;
		}

		int jjant = 1;
		for (int j = 1; j <= ff1.size(); j++)
		{
			int jj = 0;
			for (jj = jjant; jj <= ff1r.size() - 1; jj++)
			{
				if (ff1r(jj - 1) <= ff1(j - 1) && ff1r(jj) >= ff1(j - 1))
				{
					jjant = jj;
					break;
				}
			}
			ff2(j-1) = ff2r(jj-1) + (ff2r(jj) - ff2r(jj-1))*(ff1(j-1) - ff1r(jj-1)) / (ff1r(jj) - ff1r(jj-1));
		}

		fwxy[k - 1].x = ff1;
		fwxy[k - 1].y = ff2;
		
	}

	return fwxy;
}


