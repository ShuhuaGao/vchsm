#pragma once
#include <Eigen/Dense>
#include <fstream>
#include <utility>

/**
 * Serialize an Eigen matrix/vector/array 
 *@param d the Eigen variable to be serialized
 *@param os the out file stream
*/
template<typename Derived>
void serialize(const Eigen::PlainObjectBase<Derived>& d, std::ofstream& os)
{
	auto rows = d.rows(); // Eigen::Index
	auto cols = d.cols();
	//1. first write the number of rows and cols
	os.write((char*)&rows, sizeof(rows));
	os.write((char*)&cols, sizeof(cols));
	//2. write the data of the matrix/vector
	os.write((char*)d.data(), sizeof(d(0)) * rows * cols);
}


/**
 * Deserialize a file into an Eigen matrix/vector/array variable
 *@param d the Eigen variable to be deserialized into 
 *@param is the input file stream
*/
template<typename Derived>
void deserialize(Eigen::PlainObjectBase<Derived>& d, std::ifstream& is)
{
	using SizeType = decltype(std::declval<Derived>().rows());
	//1. read rows and cols
	SizeType rows{}, cols{};
	is.read(reinterpret_cast<char*>(&rows), sizeof(rows));
	is.read((char*)&cols, sizeof(cols));
	//2. read the data 
	d.resize(rows, cols);
	is.read((char*)d.data(), sizeof(d(0)) * rows * cols);
	
}