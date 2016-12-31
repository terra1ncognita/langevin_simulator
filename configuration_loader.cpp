#include "configuration_loader.h"
#include "library.h"


using json = nlohmann::json;
const double E = std::exp(1.0);
const double kBoltz = 1.38064852e-5;// (*pN um *)



/// Assign Simulation Parameters
SimulationParameters assign_simulation_parameters_from_json(SimulationParameters simp, json jsonobjsimp) {
	if (!(jsonobjsimp["expTime"].empty())) {
		simp.expTime = stod(jsonobjsimp["expTime"].get<std::string>());
	}
	if (!(jsonobjsimp["microsteps"].empty())) {
		simp.microsteps = (int)round(stod(jsonobjsimp["microsteps"].get<std::string>()));
	}
	if (!(jsonobjsimp["nTotal"].empty())) {
		simp.nTotal = (int)round(stod(jsonobjsimp["nTotal"].get<std::string>()));
	}
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
	if (!(jsonobj["ModelParameters"]["T"].empty())) {
		conf.modelParameters.T= stod(jsonobj["ModelParameters"]["T"].get<std::string>());
		conf.modelParameters.kT = kBoltz*conf.modelParameters.T;
	}
	if (!(jsonobj["ModelParameters"]["Gdepth"].empty())) {
		conf.modelParameters.G = (conf.modelParameters.kT)*stod(jsonobj["ModelParameters"]["Gdepth"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["L"].empty())) {
		conf.modelParameters.L = stod(jsonobj["ModelParameters"]["L"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["sigma"].empty())) {
		conf.modelParameters.sigma = stod(jsonobj["ModelParameters"]["sigma"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMol"].empty())) {
		conf.modelParameters.DMol = stod(jsonobj["ModelParameters"]["DMol"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DBead"].empty())) {
		conf.modelParameters.DBead = stod(jsonobj["ModelParameters"]["DBead"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMT"].empty())) {
		conf.modelParameters.DMT = stod(jsonobj["ModelParameters"]["DMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaMolAdjust"].empty())) {
		conf.modelParameters.gammaMol = conf.modelParameters.kT*stod(jsonobj["ModelParameters"]["gammaMolAdjust"].get<std::string>())/ conf.modelParameters.DMol;
	}
	if (!(jsonobj["ModelParameters"]["gammaBeadAdjust"].empty())) {
		conf.modelParameters.gammaBead = conf.modelParameters.kT*stod(jsonobj["ModelParameters"]["gammaBeadAdjust"].get<std::string>()) / conf.modelParameters.DBead;
	}
	if (!(jsonobj["ModelParameters"]["gammaMTAdjust"].empty())) {
		conf.modelParameters.gammaMT = conf.modelParameters.kT*stod(jsonobj["ModelParameters"]["gammaMTAdjust"].get<std::string>()) / conf.modelParameters.DMT;
	}
	if (!(jsonobj["ModelParameters"]["trapstiff"].empty())) {
		conf.modelParameters.trapstiff = stod(jsonobj["ModelParameters"]["trapstiff"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffL"].empty())) {
		conf.modelParameters.MTstiffL = stod(jsonobj["ModelParameters"]["MTstiffL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffR"].empty())) {
		conf.modelParameters.MTstiffR = stod(jsonobj["ModelParameters"]["MTstiffR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molstiff"].empty())) {
		conf.modelParameters.molstiff = stod(jsonobj["ModelParameters"]["molstiff"].get<std::string>());
	}
	//// Assign Initial Conditions from json
	if (!(jsonobj["InitialConditions"]["xPed"].empty())) {
		conf.initialConditions.xPed = stod(jsonobj["InitialConditions"]["xPed"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMol"].empty())) {
		conf.initialConditions.initialState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMT"].empty())) {
		conf.initialConditions.initialState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadl"].empty())) {
		conf.initialConditions.initialState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadr"].empty())) {
		conf.initialConditions.initialState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapl"].empty())) {
		conf.initialConditions.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.initialConditions.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}
	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["xMol"].empty())) {
		conf.currentState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMT"].empty())) {
		conf.currentState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadl"].empty())) {
		conf.currentState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadr"].empty())) {
		conf.currentState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
	}
	////
	return conf;
}


////Simulation params loader
SimulationParameters load_simulationparams(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json jsonsimp = fulljson["SimulationParameters"];
	SimulationParameters simp = {};
	return assign_simulation_parameters_from_json(simp, jsonsimp);
	 
}

//// Configuration creator
std::vector <Configuration> load_configuration(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json defaultjson = fulljson["Configuration"];
	Configuration default;
	default = assign_config_from_json(default, defaultjson);



	std::vector <Configuration> configurationsVector;
	Configuration iterate;

	// iterate the array of configurations
	for (json::iterator jsonit = defaultjson["Configurations"].begin(); jsonit != defaultjson["Configurations"].end(); ++jsonit) {
		iterate = default;
		iterate=assign_config_from_json(iterate, *jsonit);
		configurationsVector.push_back(iterate);
	}

	return configurationsVector;

}