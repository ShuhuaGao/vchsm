#pragma once
#include <Eigen/Dense>

namespace Eigen
{
	// Use single precision (float) as default. If necessary, change it to double. 
	using TFloat = double;

	using TMatrixX = Matrix<TFloat, -1, -1>;
	using TRowVectorX = Matrix<TFloat, 1, -1>;
	using TVectorX = Matrix<TFloat, -1, 1>; 
	using TRowVector2 = Matrix<TFloat, 1, 2>;
	using TRowVectorXc = Matrix < std::complex<TFloat>, 1, -1>;
	using TMatrixXc = Matrix<std::complex<TFloat>, -1, -1>;
	using TVectorXc = Matrix < std::complex<TFloat>, -1, 1>;
	using TRowVector4 = Matrix<TFloat, 1, 4>;;
}

enum class IndexBase
{
	Zero = 0,
	One = 1
};