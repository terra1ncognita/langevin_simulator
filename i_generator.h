#pragma once
#include <cstddef>
class IGenerator
{
public:
	virtual ~IGenerator()	{}
	virtual void generateNumbers() = 0;
	virtual const double* getNumbersBuffer() const = 0;
	virtual std::size_t getNumbersBufferSize() const = 0;
};

