#pragma once
#include <string>
# include <omp.h>
///////
char* getCmdOption(char ** begin, char ** end, const std::string & option);

///////
bool cmdOptionExists(char** begin, char** end, const std::string& option);
