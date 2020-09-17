#include "vsl_stream_wrapper.h"
#include <stdexcept>
#include <string>


VSLStreamWrapper::VSLStreamWrapper(MKL_INT generatorID, MKL_UINT seed)
{
	
	const auto creationResult = vslNewStream(&_stream, generatorID, seed);
	if (creationResult != VSL_STATUS_OK) {
		throw std::runtime_error{ std::string{ "can't create VSL stream #" } +std::to_string(generatorID) };
	}

}


VSLStreamWrapper::~VSLStreamWrapper()
{
	if (_stream != nullptr) {
		const auto deletionResult = vslDeleteStream(&_stream);
		if (deletionResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "error deleting stream" };
		}
	}
}
