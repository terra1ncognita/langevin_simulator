#pragma once
#include <string>
# include <omp.h>
#include "json.hpp"

using json = nlohmann::json;
///////
char* getCmdOption(char ** begin, char ** end, const std::string & option);

///////
bool cmdOptionExists(char** begin, char** end, const std::string& option);
////
std::string readfile(std::string filename);

json parse_json_string(std::string inputjsonstring);