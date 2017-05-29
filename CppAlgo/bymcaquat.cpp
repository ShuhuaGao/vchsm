#include "bymcaquat.h"
#include <cmath>

std::tuple<Eigen::TFloat, Eigen::TFloat, Eigen::TFloat, Eigen::TFloat> bymcaquat(Eigen::TFloat p1, Eigen::TFloat p2, Eigen::TFloat f1, Eigen::TFloat f2, Eigen::TFloat N, Eigen::TFloat fs)
{
	const Eigen::TFloat pi = 3.141592653589793;
	Eigen::TFloat w1 = 2 * pi * f1 / fs;
	Eigen::TFloat w2 = 2 * pi*f2 / fs;
	Eigen::TFloat M = round((1 / (2 * pi))*((p1 + w1*N - p2) + 0.5*N*(w2 - w1)));
	Eigen::TFloat c0 = p1;
	Eigen::TFloat c1 = w1;
	Eigen::TFloat c2 = (3 / (N*N))*(p2 - p1 - w1*N + 2 * pi*M) - (w2 - w1) / N;
	Eigen::TFloat c3 = (-2 / (N*N*N))*(p2 - p1 - w1*N + 2 * pi*M) + (w2 - w1) / (N*N);
	// since C++ only allows returning a single value, packs the four into a single tuple
	return std::make_tuple(c0, c1, c2, c3);
}
