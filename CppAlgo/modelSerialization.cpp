#include "modelSerialization.h"
#include "EigenSerialization.h"
#include <fstream>


void serializeModel(const HSMModel & model, const std::string & file)
{
	std::ofstream outFile;
	outFile.open(file, std::iostream::out | std::iostream::binary);
	auto writeStructArraySize = [&outFile](const auto& c) { // std::vector 
		auto size = c.size();
		outFile.write((char*)&size, sizeof(size));
	};

	// f0f
	serialize(model.f0f, outFile);
	// fw
	writeStructArraySize(model.fw);
	for (size_t i = 0; i < model.fw.size(); i++)
	{
		serialize(model.fw[i].x, outFile);
		serialize(model.fw[i].y, outFile);
	}
	// th and th2
	auto writeth = [&outFile, &writeStructArraySize](const auto& th) {
		writeStructArraySize(th);
		for (size_t i = 0; i < th.size(); i++)
		{
			outFile.write((char*)&th[i].a, sizeof(th[i].a)); // a scalar a
			serialize(th[i].u, outFile);
			serialize(th[i].v, outFile);
			serialize(th[i].E, outFile);
			serialize(th[i].R, outFile);
		}
	};
	// th
	writeth(model.th);
	writeth(model.th2);
	outFile.close();
}


HSMModel deserializeModel(const std::string & file)
{
	
	std::ifstream inFile;
	inFile.open(file, std::iostream::in | std::iostream::binary);
	auto readStructArraySize = [&inFile](auto& c)
	{
		decltype(c.size()) size{};
		inFile.read((char*)&size, sizeof(size));
		c.resize(size);
	};
	HSMModel model;
	// f0f
	deserialize(model.f0f, inFile);
	// fw
	readStructArraySize(model.fw);
	for (size_t i = 0; i < model.fw.size(); i++)
	{
		deserialize(model.fw[i].x, inFile);
		deserialize(model.fw[i].y, inFile);
	}
	// th and th2
	auto readth = [&inFile, &readStructArraySize](ThxyStructArray& th)
	{
		readStructArraySize(th);
		for (size_t i = 0; i < th.size(); i++)
		{
			inFile.read((char*)&th[i].a, sizeof(th[i].a));
			deserialize(th[i].u, inFile);
			deserialize(th[i].v, inFile);
			deserialize(th[i].E, inFile);
			deserialize(th[i].R, inFile);
		}
	};
	readth(model.th);
	readth(model.th2);
	return model;
}


