#pragma once
#include <algorithm>
#include <mkl_vsl.h>

class VSLStreamWrapper
{
public:
	VSLStreamWrapper(MKL_INT generatorID, MKL_UINT seed);
	// This part to prevent copying of stream from outside //////////
	VSLStreamWrapper() = default;
	VSLStreamWrapper(const VSLStreamWrapper&) = delete;
	VSLStreamWrapper& operator=(const VSLStreamWrapper&) = delete;
	VSLStreamWrapper(VSLStreamWrapper&& other) {
		std::swap(_stream, other._stream);
	}
	VSLStreamWrapper& operator=(VSLStreamWrapper&& other) {
		std::swap(_stream, other._stream);
		return *this;
	}
	//////////////////////////////////////////////////////////////////
	~VSLStreamWrapper();
	VSLStreamStatePtr get() const { return _stream; }
private:
	VSLStreamStatePtr _stream = nullptr;
};

