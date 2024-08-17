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
	const double _mean;
	const double _stDeviation;
	const std::size_t _bufferSize;
	const unsigned _threadNum;

	static std::vector <VSLStreamWrapper> _streamWrappers;
	static std::vector<double> _buffer;

	std::size_t _nPerThread;
};

