#include <memory>
#include <string>
#include <iostream>
# include <cstdlib>

#include "logger.h"



class ConsoleLogger :public ILogger
{
public:
	template <typename... Args>
	ConsoleLogger(Args&&... args) {}

	virtual void save(double xBeadl, double xMT, double xBeadr, double xMol) override {
		std::cout << xBeadl << "	" << xMT << "	" << xBeadr << "	" << xMol << endl;
	}
	virtual ~ConsoleLogger() override {

	}
};
class BinaryFileLogger :public ILogger
{
private:
	FILE* const pFile;
	int buffsize;
public:

	template <typename... Args>
	BinaryFileLogger(Args&&... args) {
		throw std::logic_error{ "wrong BinaryFileLogger construct parameter types" };
	}

	BinaryFileLogger(string filenameTemplate, string filePath)
		:pFile{ fopen((filePath + filenameTemplate + std::string{ "file.binary" }).c_str(), "wb") }
	{
		if (pFile == nullptr) {
			throw std::runtime_error{ "the file was not created" };
		}
		//std::cout << "hello" << std::endl;
	}

	virtual void save(double xBeadl, double xMT, double xBeadr, double xMol) override {
		//std::cout << xBeadl << "	" << xMT << "	" << xBeadr << "	" << xMol << endl;
		if (fwrite(&xBeadl, sizeof(double), 1, pFile) != 1) {
			throw std::runtime_error{ "not all data was written to file" };
		};
	}
	virtual ~BinaryFileLogger() override {
		fclose(pFile);
	}
};



template <typename... Args>
static std::shared_ptr <ILogger> createLogger(LoggerType type, Args&&... args) {
	if (type == LoggerType::BINARY_FILE) {
		return std::make_shared <BinaryFileLogger>(std::forward<Args>(args)...);
	}
	if (type == LoggerType::CONSOLE) {
		return std::make_shared <ConsoleLogger>(std::forward<Args>(args)...);
	}
	throw std::logic_error{ "unknown logger type" };
}