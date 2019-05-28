#include "mkl_gaussian_parallel_generator.h"
#include <stdlib.h>
#include <stdexcept>
#include <omp.h>



MklGaussianParallelGenerator::MklGaussianParallelGenerator(double mean, double stDeviation, std::size_t bufferSize, unsigned threadNum)
	:_mean{ mean }, _stDeviation{ stDeviation }, _bufferSize{bufferSize},_threadNum{threadNum}
{
	///////////////// If reproducibility from launch to launch is required seed is const, eslse seed must be random
	MKL_UINT seed = 777;
	/////////////////
	for (unsigned i = 0; i < threadNum; i++) {
		_streamWrappers.emplace_back(VSL_BRNG_MT2203 + i, seed);
	}
	_nPerThread = _bufferSize / _threadNum;
	if (_bufferSize != _nPerThread * _threadNum) {
		throw std::logic_error{ "buffsize must be multiple of number of threads" };
	}
	_buffer.resize(_bufferSize);
}



void MklGaussianParallelGenerator::generateNumbers()
{
	//// Strange that it works without shared! check this out.
#pragma omp parallel num_threads(_threadNum) default(none)
	{
		const std::size_t threadId=omp_get_thread_num();
		const std::size_t begin = _nPerThread * threadId;
		const auto generationResult = vdRngGaussian(
			VSL_RNG_METHOD_GAUSSIAN_ICDF, _streamWrappers.at(threadId).get(),
			_nPerThread, _buffer.data() + begin,
			_mean, _stDeviation
		);
		if (generationResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "can't generate numbers" };
		}
	}

	
}

const double * MklGaussianParallelGenerator::getNumbersBuffer() const
{
	return _buffer.data();
}

std::size_t MklGaussianParallelGenerator::getNumbersBufferSize() const
{
	return _bufferSize;
}
