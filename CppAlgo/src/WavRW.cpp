#include "WavRW.h"
#include "WavHeader.h"
#include <fstream>
#include <limits>

Eigen::TRowVectorX readWav(const std::string & wavFile)
{
	std::ifstream inFile;
	inFile.open(wavFile, std::ifstream::in | std::ifstream::binary);
	// read the header
	WavHeader header;
	inFile.read((char*)&header, sizeof(header));
	assert(header.num_channels == 1);
	assert(header.sample_rate == 16000);
	assert(header.bits_per_sample == 16);
	// read the data 
	auto numSample = header.subchunk2_size / 2; // 2 bytes per sample
	Eigen::Matrix<std::int16_t, 1, -1> rawData(numSample);
	inFile.read(reinterpret_cast<char*>(rawData.data()), numSample * 2);
	// normalize the data 
	inFile.close();
	constexpr auto min = std::numeric_limits<std::int16_t>::min();
	return rawData.cast<Eigen::TFloat>().array() / (-min);
}

void writeWav(Eigen::Ref<const Eigen::TRowVectorX> data, const std::string & wavFile)
{
	// first create the header
	WavHeader header;
	header.chunk_id[0] = 'R';
	header.chunk_id[1] = 'I';
	header.chunk_id[2] = 'F';
	header.chunk_id[3] = 'F';
	header.subchunk2_size = (int)data.size() * 2;
	header.chunk_size = 36 + header.subchunk2_size;
	header.format[0] = 'W';
	header.format[1] = 'A';
	header.format[2] = 'V';
	header.format[3] = 'E';
	header.subchunk1_id[0] = 'f';
	header.subchunk1_id[1] = 'm';
	header.subchunk1_id[2] = 't';
	header.subchunk1_id[3] = 32;
	header.subchunk1_size = 16;
	header.audio_format = 1;
	header.num_channels = 1;
	header.sample_rate = 16000;
	header.byte_rate = 32000;
	header.block_align = 2;
	header.bits_per_sample = 16;
	header.subchunk2_id[0] = 'd';
	header.subchunk2_id[1] = 'a';
	header.subchunk2_id[2] = 't';
	header.subchunk2_id[3] = 'a';

	// write the header
	std::ofstream outFile;
	outFile.open(wavFile, std::iostream::out | std::iostream::binary);
	outFile.write((char*)&header, sizeof(header));

	// transform the normalized data to 16bits
	constexpr auto min = std::numeric_limits<std::int16_t>::min();
	auto y = (-min * data).cast<std::int16_t>().eval();
	// write the data
	outFile.write(reinterpret_cast<char*>(y.data()), header.subchunk2_size);
	outFile.close();
}
