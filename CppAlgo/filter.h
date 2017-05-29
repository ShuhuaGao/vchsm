#pragma once
#include <Eigen/Dense>

/**
 * 1-D digital filter
 *@param[out] speech the data to be filtered in a vector form
 *@remark This is similar but limited implementation of MATLAB filter function, see https://www.mathworks.com/help/matlab/ref/filter.html.
 *@remark Only for the specific usage in the mfcc.m: speech = filter( [1 -alpha], 1, speech )
 *@note This is a rational transfer function based filter, whose equation is y(n) = x(n) - alpha*x(n-1).
*/
template<typename Derived>
void filter(Eigen::TFloat alpha, Eigen::MatrixBase<Derived>& speech)
{
	auto size = speech.size();
	assert(size > 0);
	speech.tail(size - 1).reverse() -= alpha * speech.head(size - 1).reverse();
}

/**
* 1-D digital filter
*@param[out] ye the data to be filtered in a vector form
*@remark This is similar a but limited implementation of MATLAB filter function, see https://www.mathworks.com/help/matlab/ref/filter.html.
*@remark Only for the specific usage in (1) stochasticanalysis.m ye=filter(fir1(48,500/(0.5*fs),'high'),1,ye);  and 
*(2) HSManalyze.m ye=filter(fir1(48,500/(0.5*fs),'high'),1,ye); 
*@note Since fs is fixed to be 16000, we can precompute fir1(48,500/(0.5*fs),'high') and then finds the implementation of filter here. 
*/
template<typename Derived>
void filter(Eigen::MatrixBase<Derived>& ye)
{
	assert(ye.size() >= 49);
	// fir1(48,500/(0.5*fs),'high')
	int nb = 48;
	Eigen::TRowVectorX b(49);
	b << 0.001059931368409, 0.001138125306749, 0.001277574223688, 0.001448045264934, 0.001592223730547, 0.001627664902621, 0.001451645795686, 0.000948655673462, -0.000000000000000, -0.001505233343453, -0.003658762031140, -0.006522255494314, -0.010118042885623, -0.014422295227610, -0.019361418183423, -0.024812144650090, -0.030605518262808, -0.036534635299900, -0.042365691542421, -0.047851590755281, -0.052747138661000, -0.056824691882445, -0.059889069935446, -0.061790576460614, 0.936526668835176, -0.061790576460614, -0.059889069935446, -0.056824691882445, -0.052747138661000, -0.047851590755281, -0.042365691542421, -0.036534635299900, -0.030605518262808, -0.024812144650090, -0.019361418183423, -0.014422295227610, -0.010118042885623, -0.006522255494314, -0.003658762031140, -0.001505233343453, -0.000000000000000, 0.000948655673462, 0.001451645795686, 0.001627664902621, 0.001592223730547, 0.001448045264934, 0.001277574223688, 0.001138125306749, 0.001059931368409;
	// since a = 1 in the above filter function, the algorithm is: y(n) = b(1)x(n) + b(2)x(n-1) +... + b(nb+1)x(n-nb)
	// since we filter in place, we have to in a backward way 
	for (int n = (int)ye.size() - 1; n >= nb; n--)
	{
		ye(n) = b.dot(ye.segment(n - nb, nb + 1).reverse());
	}
	// when n < nb, there is not enough source data, and use 0 
	for (int n = nb - 1; n >= 0; n--)
	{
		ye(n) = b.head(n + 1).dot(ye.head(n + 1).reverse());
	}

}



