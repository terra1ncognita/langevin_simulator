#include "configuration_loader.h"
#include "library.h"
#include<fstream>

using json = nlohmann::json;
const double kBoltz = 1.38064852e-5;// (*pN um *)

// Assign Simulation Parameters
SimulationParameters assign_simulation_parameters_from_json(SimulationParameters simp, json jsonobjsimp) {
	if (!(jsonobjsimp["expTime"].empty())) {
		simp.expTime = stod(jsonobjsimp["expTime"].get<std::string>());
	}
	if (!(jsonobjsimp["simulationTime"].empty())) {
		simp.simulationTime = stod(jsonobjsimp["simulationTime"].get<std::string>());
	}
	if (!(jsonobjsimp["iterationsbetweenSavings"].empty())) {
		simp.iterationsbetweenSavings = stoi(jsonobjsimp["iterationsbetweenSavings"].get<std::string>());
	}
	if (!(jsonobjsimp["iterationsbetweenTrapsUpdate"].empty())) {
		simp.iterationsbetweenTrapsUpdate = stoi(jsonobjsimp["iterationsbetweenTrapsUpdate"].get<std::string>());
	}

	if (!(jsonobjsimp["totalsavings"].empty())) {
		simp.totalsavings = stoi(jsonobjsimp["totalsavings"].get<std::string>());
	}
	else {
		simp.totalsavings = int((simp.simulationTime / simp.expTime) / simp.iterationsbetweenSavings);
	}

	if (!(jsonobjsimp["buffsize"].empty())) {
		simp.buffsize = stoi(jsonobjsimp["buffsize"].get<std::string>());
	}
	if (!(jsonobjsimp["randomsPeriter"].empty())) {
		simp.randomsPeriter = stoi(jsonobjsimp["randomsPeriter"].get<std::string>());
	}

	if (!(jsonobjsimp["stepsperbuffer"].empty())) {
		simp.stepsperbuffer = stoi(jsonobjsimp["stepsperbuffer"].get<std::string>());
	}
	else {
		simp.stepsperbuffer = static_cast<int>(std::floor(simp.buffsize / simp.randomsPeriter));
	}

	if (simp.stepsperbuffer % simp.iterationsbetweenSavings != 0) {
		throw std::runtime_error{ "Please check that totalsavings/stepsperbuffer is integer" };
	}
	simp.savingsPerMacrostep = simp.stepsperbuffer / simp.iterationsbetweenSavings;

	if (simp.totalsavings % simp.savingsPerMacrostep != 0) {
		throw std::runtime_error{ "Please check that totalsavings/stepsperbuffer is integer" };
	}
	simp.macrostepMax = simp.totalsavings / simp.savingsPerMacrostep;

	if (simp.iterationsbetweenTrapsUpdate % simp.iterationsbetweenSavings != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenTrapsUpdate/iterationsbetweenSavings is integer" };
	}
	simp.trapsUpdateTest = simp.iterationsbetweenTrapsUpdate / simp.iterationsbetweenSavings;

	return simp;
}

// Assign configuration from json object
Configuration assign_config_from_json(Configuration conf, json jsonobj) {
	//// Assign Logger Parameters from json
	if (!(jsonobj["Name"].empty())) {
		conf.loggerParameters.name = jsonobj["Name"].get<std::string>();
	}
	if (!(jsonobj["LoggerParameters"]["FilePath"].empty())) {
		conf.loggerParameters.filepath= jsonobj["LoggerParameters"]["FilePath"].get<std::string>();
	}
	
	//// Assign Model Parameters from json
	
	if (!(jsonobj["Name"].empty())) {
		conf.modelParameters.name = jsonobj["Name"].get<std::string>();
	}

	if (!(jsonobj["ModelParameters"]["T"].empty())) {
		conf.modelParameters.T= stod(jsonobj["ModelParameters"]["T"].get<std::string>());
		conf.modelParameters.kT = kBoltz*conf.modelParameters.T;
	}
	if (!(jsonobj["ModelParameters"]["Gdepth"].empty())) {
		conf.modelParameters.G = (conf.modelParameters.kT)*stod(jsonobj["ModelParameters"]["Gdepth"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["Gdepth2"].empty())) {
		conf.modelParameters.G2 = (conf.modelParameters.kT)*stod(jsonobj["ModelParameters"]["Gdepth2"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["L"].empty())) {
		conf.modelParameters.L = stod(jsonobj["ModelParameters"]["L"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["sigma"].empty())) {
		conf.modelParameters.sigma = stod(jsonobj["ModelParameters"]["sigma"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["A"].empty())) {
		conf.modelParameters.A = stod(jsonobj["ModelParameters"]["A"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["m"].empty())) {
		conf.modelParameters.m = stod(jsonobj["ModelParameters"]["m"].get<std::string>());
	}


	if (!(jsonobj["ModelParameters"]["gammaMol"].empty())) {
		conf.modelParameters.gammaMol = stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
		conf.modelParameters.DMol = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
	}
	

	// Molecular stiffness
	if (!(jsonobj["ModelParameters"]["molStiffWeakSlope"].empty())) {
		conf.modelParameters.molStiffWeakSlope = stod(jsonobj["ModelParameters"]["molStiffWeakSlope"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molStiffBoundary"].empty())) {
		conf.modelParameters.molStiffBoundary = stod(jsonobj["ModelParameters"]["molStiffBoundary"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molStiffStrongSlope"].empty())) {
		conf.modelParameters.molStiffStrongSlope = stod(jsonobj["ModelParameters"]["molStiffStrongSlope"].get<std::string>());
	}
	
	
	if (!(jsonobj["ModelParameters"]["kOn1"].empty())) {
		conf.modelParameters.kOn1 = stod(jsonobj["ModelParameters"]["kOn1"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kOff1"].empty())) {
		conf.modelParameters.kOff1 = stod(jsonobj["ModelParameters"]["kOff1"].get<std::string>());
	}
	/*if (!(jsonobj["ModelParameters"]["kOn2"].empty())) {
		conf.modelParameters.kOn2 = stod(jsonobj["ModelParameters"]["kOn2"].get<std::string>());
	}*/
	/*if (!(jsonobj["ModelParameters"]["kOff2"].empty())) {
		conf.modelParameters.kOff2 = stod(jsonobj["ModelParameters"]["kOff2"].get<std::string>());
	}*/

	if (!(jsonobj["ModelParameters"]["rotFriction"].empty())) {
		conf.modelParameters.rotFriction = stod(jsonobj["ModelParameters"]["rotFriction"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["rotStiffness"].empty())) {
		conf.modelParameters.rotStiffness = stod(jsonobj["ModelParameters"]["rotStiffness"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molLength"].empty())) {
		conf.modelParameters.molLength = stod(jsonobj["ModelParameters"]["molLength"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["domainsDistance"].empty())) {
		conf.modelParameters.domainsDistance = stod(jsonobj["ModelParameters"]["domainsDistance"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["rotWellWidth"].empty())) {
		conf.modelParameters.rotWellWidth = stod(jsonobj["ModelParameters"]["rotWellWidth"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["rotWellDepth"].empty())) {
		conf.modelParameters.rotWellDepth = (conf.modelParameters.kT) * stod(jsonobj["ModelParameters"]["rotWellDepth"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["iniPhi"].empty())) {
		conf.modelParameters.iniPhi = stod(jsonobj["ModelParameters"]["iniPhi"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["bindingDynamics"].empty())) {
		conf.modelParameters.bindingDynamics = bool(stoi(jsonobj["ModelParameters"]["bindingDynamics"].get<std::string>()));
	}

	if (!(jsonobj["ModelParameters"]["B"].empty())) {
		conf.modelParameters.B = stod(jsonobj["ModelParameters"]["B"].get<std::string>());
	}


	//// Assign Initial Conditions from json
	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["xMol"].empty())) {
		conf.initialConditions.initialState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
		conf.currentState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["direction"].empty())) {
		conf.currentState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
		conf.initialConditions.initialState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["phi"].empty())) {
		conf.currentState.phi = stod(jsonobj["InitialConditions"]["phi"].get<std::string>());
		conf.initialConditions.initialState.phi = stod(jsonobj["InitialConditions"]["phi"].get<std::string>());
	}

	
	//TODO: generalize to arbitrary Markov chain
	if (conf.modelParameters.numStates != 2) {
		throw std::runtime_error{ "Current implementation is only for 2 states!" };
	}

	conf.modelParameters.transitionMatrix = new double*[conf.modelParameters.numStates];
	for (int i = 0; i < conf.modelParameters.numStates; ++i) {
		conf.modelParameters.transitionMatrix[i] = new double[conf.modelParameters.numStates];
	}

	conf.modelParameters.kOff1 = 0.0;

	conf.modelParameters.transitionMatrix[0][0] = -conf.modelParameters.kOn1;
	conf.modelParameters.transitionMatrix[0][1] = conf.modelParameters.kOn1;

	conf.modelParameters.transitionMatrix[1][1] = conf.modelParameters.kOff1;
	conf.modelParameters.transitionMatrix[1][2] = -conf.modelParameters.kOff1;

	return conf;
}


//Simulation params loader
SimulationParameters load_simulationparams(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json jsonsimp = fulljson["SimulationParameters"];
	SimulationParameters simp;
	return assign_simulation_parameters_from_json(simp, jsonsimp);
}

//// Configuration creator new
std::vector <Configuration> load_configuration(std::string paramInputFilename, unsigned nThreads) {

	std::ifstream fin(paramInputFilename);
	std::vector<std::string> configs;
	std::copy(std::istream_iterator<std::string>(fin),
		std::istream_iterator<std::string>(),
		std::back_inserter(configs));
	fin.close();

	//std::vector <std::string> configs = split(paramInputFilename, ";");
	if (configs.size() % nThreads != 0) {
		throw std::runtime_error{ "Number of configs is not multiple of number of cores" };
	}
	std::vector <Configuration> configurationsVector;
	for (const auto& config : configs) {
		json fulljson = parse_json_string(readfile(config));
		json defaultjson = fulljson["Configuration"];
		Configuration iterate;
		iterate = assign_config_from_json(iterate, defaultjson);
		configurationsVector.push_back(iterate);
	}
	return configurationsVector;
}
