#include "vec2frames.h"
#include "concat.h"
#include "seq.h"
#include "hamming.h"



Eigen::TMatrixX vec2frames(Eigen::Ref<const Eigen::TVectorX>vec, int Nw, int Ns)
{
	auto L = vec.size();
	auto M = static_cast<int>(std::floor((L - Nw) / Eigen::TFloat(Ns) + 1));

	// case 'cols'
	auto indf = seq<Eigen::RowVectorXi>(0, M - 1) * Ns;
	auto inds = seq<Eigen::VectorXi>(1, Nw);
	// indexes = indf(ones(Nw,1),:) + inds(:,ones(1,M)); MATLAB
	auto indexes = indf.replicate(Nw, 1) + inds.replicate(1, M);

	// MATALB: frames = vec( indexes ); use the above indexes matrix to index the vec into a same-dimensional matrix
	// NO direct support in Eigen. However, here notice that each column of indexes is a sequential sequence and vel is a column vector.
	auto rowNum = indexes.rows();
	auto colNum = indexes.cols();
	Eigen::TMatrixX frames(rowNum, colNum);
	for (Eigen::Index i = 0; i < frames.cols(); i++)
	{
		frames.col(i) = vec.segment(indexes(0, i) - 1, rowNum); // Eigen starts indexing from 0
	}

	// generate window samples with hamming
	auto window = hamming(Nw);

	if (window.size() == Nw)
	{
		frames = window.asDiagonal() * frames;
	}
	return frames;
}
