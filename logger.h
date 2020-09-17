#pragma once
#include <memory>

// Logger of simulation data
class ILogger
{
public:
	virtual ~ILogger() {}
	virtual void save(double xBeadl, double xMT, double xBeadr, double xMol) = 0;
};


enum class LoggerType {
	CONSOLE,
	BINARY_FILE
};

template <typename... Args>
static std::shared_ptr <ILogger> createLogger(LoggerType type, Args&&... args);