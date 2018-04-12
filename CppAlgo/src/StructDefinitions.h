#pragma once
#include "types.h"
#include "PicosElement.h"
#include <vector>
/** \file StuctDefinitions.h
 * \brief In MATLAB functions, there are some struct definitions. We list these definitions here.
*/

//!
//! \struct Structure to be used in variable aux in clustkmeans.m
struct AuxElementType
{
	int act = 0;
	Eigen::RowVectorXi ind;
	Eigen::TRowVectorX sumd;
};

/**
 * \struct Structure to be used for variable thxy in HSMptraining.m
*/
struct ThxyElementType
{
	Eigen::TFloat a = 0;
	Eigen::TVectorX u;
	Eigen::TVectorX v;
	Eigen::TMatrixX E;
	Eigen::TMatrixX R;

};
using ThxyStructArray = std::vector<ThxyElementType>;

/**
* \struct Structure to be used for variable fwxy and fwyx in automaticth2wfw.m
*/
struct FwxyElementType
{
	Eigen::RowVectorXi x;
	Eigen::TRowVectorX y;

};
using FwxyStructArray = std::vector<FwxyElementType>;

struct PolosElementType
{
	Eigen::TRowVectorX ai1, ai2, cc1, cc2;
	Eigen::TRowVectorXc p1, p2;
};

using PolosStructArray = std::vector<PolosElementType>;


/**
*\struct Structure for the final model parameters.
*/
struct HSMModel
{
	ThxyStructArray th;
	ThxyStructArray th2;
	FwxyStructArray fw;
	Eigen::TRowVector4 f0f;

	// default constructor
	HSMModel()
	{}

	// copy constructor 
	HSMModel(const HSMModel& other)
		: th{ other.th }, th2{ other.th2 }, fw{ other.fw }, f0f{ other.f0f }
	{}

	// move constructor
	HSMModel(HSMModel&& other) // Eigen objects have no move constructor yet.
		: th{ std::move(other.th) }, th2{ std::move(other.th2) }, fw{ std::move(other.fw) }, f0f{ std::move(other.f0f) }
	{}

	// copy assignment
	HSMModel& operator=(const HSMModel& other)
	{
		th = other.th;
		th2 = other.th2;
		fw = other.fw;
		f0f = other.f0f;
		return *this;
	}
	// move assignment
	HSMModel& operator=(HSMModel&& other)
	{
		th = std::move(other.th);
		th2 = std::move(other.th2);
		fw = std::move(other.fw);
		f0f = std::move(other.f0f);
		return *this;
	}
};

using PicosStructArray = std::vector<PicosElement>;


/**
 * \struct Structure for the thz variable in HSMptraining.m
*/
struct ThzElementType
{
	Eigen::TFloat a = 0;
	Eigen::TVectorX u;
	Eigen::TMatrixX E;
};

using ThzStructArray = std::vector<ThzElementType>;




