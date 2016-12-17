#pragma once
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "json.hpp"
# include <omp.h>
#include "configuration.h"

using json = nlohmann::json;

std::string readfile(std::string filename);

json parse_json_string(std::string inputjsonstring);
SimulationParameters assign_simulation_parameters_from_json(SimulationParameters simp, json jsonobjsimp);

Configuration assign_config_from_json(Configuration conf, json jsonobj);

SimulationParameters load_simulationparams(std::string paramInputFilename);
std::vector <Configuration> load_configuration(std::string paramInputFilename);