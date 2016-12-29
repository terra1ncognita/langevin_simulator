#pragma once
#include "i_generator.h"
#include "vsl_stream_wrapper.h"
#include <vector>
class MklGaussianParallelGenerator :
	public IGenerator
{
public:
	MklGaussianParallelGenerator(double mean, double stDeviation, std::size_t bufferSize, unsigned threadNum);
	virtual void generateNumbers() override;
	virtual const double* getNumbersBuffer() const override;
	virtual std::size_t getNumbersBufferSize() const override;
private:
	double _mean;
	double _stDeviation;
	std::size_t _bufferSize;
	unsigned _threadNum;
	std::vector <VSLStreamWrapper> _streamWrappers;
	std::vector<double> _buffer;
	std::size_t _nPerThread;
	

};

