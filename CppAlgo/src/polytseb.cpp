#include "polytseb.h"

std::pair<TMatrixX, TMatrixX> polytseb(int K)
{
	TMatrixX T = TMatrixX::Zero(K + 1, K + 1);
	// Note the index in Eigen (C++) starts from 0 while MATLAB begins from 1
	T(0, K) = 1;
	T(1, K - 1) = 1;
	TMatrixX U = TMatrixX::Zero(K + 1, K + 1);
	U(0, K) = 1;
	U(1, K - 1) = 2;

	for (int k = 2; k <= K; k++)
	{
		// Note again: in C++ the index starts from 0
		// T(k + 1, :) = 2 * [T(k, 2:K + 1) 0] - T(k - 1, :); (MATLAB)
		// Hint: T(k, 2:K + 1) gets the 2 to K + 1 elements of the kth row
		T.row(k) = 2 * (TRowVectorX(K + 1) << T.row(k - 1).segment(1, K), 0).finished() - T.row(k - 2);
		// U(k+1,:)=[U(k,2:K+1) 0]+T(k+1,:);  (MATLAB)
		U.row(k) = (TRowVectorX(K + 1) << U.row(k - 1).segment(1, K), 0).finished() + T.row(k);
	}
	// T=T(2:K+1,:); U=U(1:K,:); (MATLAB)
	// Note: the bottom-most K rows of T and the top-most K rows of U are extracted
	// Since C++ only return a single value, we have to pack T and U into a pair.
	return std::make_pair(T.bottomRows(K), U.topRows(K));
}





