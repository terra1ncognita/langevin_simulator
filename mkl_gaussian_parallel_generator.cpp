#include "mkl_gaussian_parallel_generator.h"
#include <stdlib.h>
#include <stdexcept>
#include <omp.h>
#include <ctime>
#include <iostream>


MklGaussianParallelGenerator::MklGaussianParallelGenerator(double mean, double stDeviation, std::size_t bufferSize, unsigned threadNum)
	:_mean{ mean }, _stDeviation{ stDeviation }, _bufferSize{bufferSize}, _threadNum{threadNum}
{
	int factor = 10;

	_nPerThread = _bufferSize / (factor * _threadNum );
	if (_bufferSize != _nPerThread * _threadNum * factor) {
		throw std::logic_error{ "buffsize must be multiple of number of threads" };
	}
	std::cout << "Random numbers buffer alocation of " << _bufferSize << " bytes started...  ";
	_buffer.resize(_bufferSize);
	std::cout << "Done" << std::endl;

	///////////////// If reproducibility from launch to launch is required, then seed is const, else seed must be random

	// MKL_UINT seed = 777;
	srand(time(NULL));

	/////////////////

	MKL_UINT seed = static_cast<MKL_UINT>(rand());
	_sw = VSLStreamWrapper(VSL_BRNG_MT2203, seed);
	std::cout << "Added generator thread with seed " << seed << std::endl;

	for (unsigned i = 0; i < threadNum; i++) {
		MKL_UINT seed = static_cast<MKL_UINT>(rand());
		_streamWrappers.emplace_back(VSL_BRNG_MT2203 + i, seed);
		std::cout << "Added generator thread with seed " << seed << std::endl;
	}
}

void MklGaussianParallelGenerator::generateNumbers()
{
	#pragma omp parallel for num_threads(_threadNum) default(none) shared(_sw) schedule(dynamic)
	for (int i = 0; i < _bufferSize; i += _nPerThread) {
		const auto begin = _buffer.data() + i;
		const auto generationResult = vdRngGaussian(
			VSL_RNG_METHOD_GAUSSIAN_ICDF, _streamWrappers.at(omp_get_thread_num()).get(),
			_nPerThread, begin, _mean, _stDeviation
		);
		if (generationResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "can't generate numbers" };
		}
	}

	//std::cout << _buffer[0] << " " << _buffer[1] << " " << _buffer[_nPerThread] << " " << _buffer[_nPerThread + 1] <<  std::endl;
	//// Strange that it works without shared! check this out.
	//#pragma omp parallel num_threads(_threadNum) default(none)
	//{
	//	const std::size_t threadId = omp_get_thread_num();
	//	const auto begin = _buffer.data() + _nPerThread * threadId;
	//	const auto generationResult = vdRngGaussian(
	//		VSL_RNG_METHOD_GAUSSIAN_ICDF, _streamWrappers.at(threadId).get(),
	//		_nPerThread, begin, _mean, _stDeviation
	//	);
	//	if (generationResult != VSL_STATUS_OK) {
	//		throw std::runtime_error{ "can't generate numbers" };
	//	}
	//}	
}

const double * MklGaussianParallelGenerator::getNumbersBuffer() const
{
	return _buffer.data();
}

std::size_t MklGaussianParallelGenerator::getNumbersBufferSize() const
{
	return _bufferSize;
}
