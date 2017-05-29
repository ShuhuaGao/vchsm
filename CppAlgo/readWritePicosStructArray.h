#pragma once
#include "mex_headers.h"
#include "PicosElement.h"
#include "MATLAB_Eigen.h"
#include "StructDefinitions.h"
#include <vector>

/**
 * Read a structure array (picos) into a C++ vector of structures
 *@param pArray MATLAB structure array
 *@return a C++ vector of structures
*/
inline PicosStructArray readPicosFromMatlab(const mxArray* pArray)
{
	assert(mxIsStruct(pArray));
	auto n = mxGetNumberOfElements(pArray); // length of struct elements
	PicosStructArray picos(n);
	for (int i = 0; i < n; i++)
	{
		picos[i].e = MatlabToEigen<Eigen::RowVectorXd>(mxGetField(pArray, i, "e")).cast<Eigen::TFloat>();
		picos[i].pm = MatlabToScalar<int>(mxGetField(pArray, i, "pm"));
		picos[i].a = MatlabToEigen<Eigen::RowVectorXd>(mxGetField(pArray, i, "a")).cast<Eigen::TFloat>();
		picos[i].alfa = MatlabToScalar<Eigen::TFloat>(mxGetField(pArray, i, "alfa"));
		
		picos[i].f0 = MatlabToScalar<Eigen::TFloat>(mxGetField(pArray, i, "f0"));
		picos[i].mfcc = MatlabToEigen<Eigen::RowVectorXd>(mxGetField(pArray, i, "mfcc")).cast<Eigen::TFloat>();
		picos[i].p = MatlabToEigen<Eigen::RowVectorXd>(mxGetField(pArray, i, "p")).cast<Eigen::TFloat>();
		
	}
	return picos;
}


/**
 * Write picos as a struct array to MATLAB
 *@param picos the original data 
 *@return the struct array created for MATLAB
*/
inline mxArray* writePicosToMatlab(const PicosStructArray& picos)
{
	auto n = picos.size();
	const char* fieldNames[] = { "a", "alfa", "e", "f0", "mfcc", "p", "pm" };
	// create a structure array 
	auto pPicosArray = mxCreateStructMatrix(1, n, 7, fieldNames);
	// set the field for each element
	for (int i = 0; i < n; i++)
	{
		// first create the field and then set it to an element
		mxSetField(pPicosArray, i, "a", EigenToMatlab(picos[i].a));
		mxSetField(pPicosArray, i, "alfa", ScalarToMatlab(picos[i].alfa));
		mxSetField(pPicosArray, i, "e", EigenToMatlab(picos[i].e));
		mxSetField(pPicosArray, i, "f0", ScalarToMatlab(picos[i].f0));
		mxSetField(pPicosArray, i, "mfcc", EigenToMatlab(picos[i].mfcc));
		mxSetField(pPicosArray, i, "p", EigenToMatlab(picos[i].p));
		mxSetField(pPicosArray, i, "pm", ScalarToMatlab(picos[i].pm));
	}

	return pPicosArray;
}
