#pragma once
#include "types.h"
#include <vector>

class RealTimeConverter
{
public:
	RealTimeConverter();
	~RealTimeConverter();
private:
	int mFragmentSize;
	int mOverlapSize;
	std::vector<Eigen::TFloat> buffer;

public:
	int fragmentSize() const
	{
		return mFragmentSize;
	}

	void setFragmentSize(int fragmentSize)
	{
		mFragmentSize = fragmentSize;
	}

	int overlapSize() const
	{
		return mOverlapSize;
	}

	void setOverlapSize(int overlapSize)
	{
		mOverlapSize = overlapSize;
	}

public:
	void convertNewFragment(const Eigen::Ref<Eigen::TRowVectorX> x);
	Eigen::TRowVectorX getAssembledConversion(); 
};

