#include "library.h"
#include <algorithm>
#include <fstream>

// Functions to read from the command line string parameters
///////
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return 0;
}
///////
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

///////////

//// Assisting functions for configuration loading
//function to read files
std::string readfile(std::string filename) {
	std::string filecontents = "";
	std::string line;
	std::ifstream myfile(filename);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			filecontents = filecontents + line;
		}
		myfile.close();
	}

	else std::cout << "Unable to open file";
	return filecontents;
}
// function to parse json objects
json parse_json_string(std::string inputjsonstring) {
	std::stringstream ss;
	ss << inputjsonstring;
	json j = json::parse(ss);
	//another way to do it is to use inputjsonstring.c_str()
	return j;

}
// function to split string
std::vector<std::string> split(const std::string& str, const std::string& delim)
{
	std::vector<std::string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == std::string::npos){
			pos = str.length();
		}
		std::string token = str.substr(prev, pos - prev);
		if (!token.empty()) {
			tokens.push_back(token);
		}
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}