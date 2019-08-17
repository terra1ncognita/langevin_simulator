#include "configuration_loader.h"
#include "library.h"


using json = nlohmann::json;
const double kBoltz = 1.38064852e-5;// (*pN um *)

/// Assign Simulation Parameters
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
		simp.totalsavings = int(simp.simulationTime / simp.iterationsbetweenSavings);
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

	if (simp.iterationsbetweenSavings % simp.stepsperbuffer != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenSavings/stepsperbuffer is integer" };
	}
	if (simp.iterationsbetweenTrapsUpdate % simp.iterationsbetweenSavings != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenTrapsUpdate/iterationsbetweenSavings is integer" };
	}

	simp.macrostepMax = simp.iterationsbetweenSavings / simp.stepsperbuffer;
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
	if (!(jsonobj["ModelParameters"]["gammaBeadL"].empty())) {
		conf.modelParameters.gammaBeadL = stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
		conf.modelParameters.DBeadL = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaBeadR"].empty())) {
		conf.modelParameters.gammaBeadR = stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
		conf.modelParameters.DBeadR = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaMT"].empty())) {
		conf.modelParameters.gammaMT =stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
		conf.modelParameters.DMT = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMol"].empty())) {
		conf.modelParameters.DMol = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DBeadL"].empty())) {
		conf.modelParameters.DBeadL = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DBeadR"].empty())) {
		conf.modelParameters.DBeadR = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMT"].empty())) {
		conf.modelParameters.DMT = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["trapstiffL"].empty())) {
		conf.modelParameters.trapstiffL = stod(jsonobj["ModelParameters"]["trapstiffL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["trapstiffR"].empty())) {
		conf.modelParameters.trapstiffR = stod(jsonobj["ModelParameters"]["trapstiffR"].get<std::string>());
	}

	// Friction coeffitient for quasiviscous molecule - MT interaction
	if (!(jsonobj["ModelParameters"]["gammaQuasiviscous"].empty())) {
		conf.modelParameters.gammaQuasiviscous = stod(jsonobj["ModelParameters"]["gammaQuasiviscous"].get<std::string>());
	}

	// MT stiffness Left
	if (!(jsonobj["ModelParameters"]["MTstiffWeakSlopeL"].empty())) {
		conf.modelParameters.MTstiffWeakSlopeL = stod(jsonobj["ModelParameters"]["MTstiffWeakSlopeL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffWeakBoundaryL"].empty())) {
		conf.modelParameters.MTstiffWeakBoundaryL = stod(jsonobj["ModelParameters"]["MTstiffWeakBoundaryL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicAL"].empty())) {
		conf.modelParameters.MTstiffParabolicAL = stod(jsonobj["ModelParameters"]["MTstiffParabolicAL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicBL"].empty())) {
		conf.modelParameters.MTstiffParabolicBL = stod(jsonobj["ModelParameters"]["MTstiffParabolicBL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicCL"].empty())) {
		conf.modelParameters.MTstiffParabolicCL = stod(jsonobj["ModelParameters"]["MTstiffParabolicCL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongBoundaryL"].empty())) {
		conf.modelParameters.MTstiffStrongBoundaryL = stod(jsonobj["ModelParameters"]["MTstiffStrongBoundaryL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongSlopeL"].empty())) {
		conf.modelParameters.MTstiffStrongSlopeL = stod(jsonobj["ModelParameters"]["MTstiffStrongSlopeL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongIntersectL"].empty())) {
		conf.modelParameters.MTstiffStrongIntersectL = stod(jsonobj["ModelParameters"]["MTstiffStrongIntersectL"].get<std::string>());
	}

	// MT stiffness Right
	if (!(jsonobj["ModelParameters"]["MTstiffWeakSlopeR"].empty())) {
		conf.modelParameters.MTstiffWeakSlopeR = stod(jsonobj["ModelParameters"]["MTstiffWeakSlopeR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffWeakBoundaryR"].empty())) {
		conf.modelParameters.MTstiffWeakBoundaryR = stod(jsonobj["ModelParameters"]["MTstiffWeakBoundaryR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicAR"].empty())) {
		conf.modelParameters.MTstiffParabolicAR = stod(jsonobj["ModelParameters"]["MTstiffParabolicAR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicBR"].empty())) {
		conf.modelParameters.MTstiffParabolicBR = stod(jsonobj["ModelParameters"]["MTstiffParabolicBR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicCR"].empty())) {
		conf.modelParameters.MTstiffParabolicCR = stod(jsonobj["ModelParameters"]["MTstiffParabolicCR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongBoundaryR"].empty())) {
		conf.modelParameters.MTstiffStrongBoundaryR = stod(jsonobj["ModelParameters"]["MTstiffStrongBoundaryR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongSlopeR"].empty())) {
		conf.modelParameters.MTstiffStrongSlopeR = stod(jsonobj["ModelParameters"]["MTstiffStrongSlopeR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongIntersectR"].empty())) {
		conf.modelParameters.MTstiffStrongIntersectR = stod(jsonobj["ModelParameters"]["MTstiffStrongIntersectR"].get<std::string>());
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
	
	
	if (!(jsonobj["ModelParameters"]["MTlength"].empty())) {
		conf.modelParameters.MTlength = stod(jsonobj["ModelParameters"]["MTlength"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molstiff"].empty())) {
		conf.modelParameters.molstiff = stod(jsonobj["ModelParameters"]["molstiff"].get<std::string>());
	}

	if (!(jsonobj["ModelParameters"]["feedbackFreq"].empty())) {
		conf.modelParameters.feedbackFreq = stod(jsonobj["ModelParameters"]["feedbackFreq"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DmblMoveAmplitude"].empty())) {
		conf.modelParameters.DmblMoveAmplitude = stod(jsonobj["ModelParameters"]["DmblMoveAmplitude"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["prestretchTotalForce"].empty())) {
		conf.modelParameters.prestretchTotalForce = stod(jsonobj["ModelParameters"]["prestretchTotalForce"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["movementTotalForce"].empty())) {
		conf.modelParameters.movementTotalForce = stod(jsonobj["ModelParameters"]["movementTotalForce"].get<std::string>());
	}

	if (!(jsonobj["ModelParameters"]["kOn1"].empty())) {
		conf.modelParameters.kOn1 = stod(jsonobj["ModelParameters"]["kOn1"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kOff1"].empty())) {
		conf.modelParameters.kOff1 = stod(jsonobj["ModelParameters"]["kOff1"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kOn2"].empty())) {
		conf.modelParameters.kOn2 = stod(jsonobj["ModelParameters"]["kOn2"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kOff2"].empty())) {
		conf.modelParameters.kOff2 = stod(jsonobj["ModelParameters"]["kOff2"].get<std::string>());
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
		conf.initialConditions.initialState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.initialConditions.initialState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["direction"].empty())) {
		conf.currentState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
		conf.initialConditions.initialState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
	}
	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["xTrapl"].empty())) {
		conf.currentState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.currentState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}

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

	//TODO: generalize to arbitrary Markov chain
	if (conf.modelParameters.numStates != 3) {
		throw std::runtime_error{ "Current implementation is only for 3 states!" };
	}

	conf.modelParameters.transitionMatrix = new double*[conf.modelParameters.numStates];
	for (int i = 0; i < conf.modelParameters.numStates; ++i) {
		conf.modelParameters.transitionMatrix[i] = new double[conf.modelParameters.numStates];
	}

	conf.modelParameters.transitionMatrix[0][0] = -conf.modelParameters.kOn1;
	conf.modelParameters.transitionMatrix[0][1] = conf.modelParameters.kOn1;
	conf.modelParameters.transitionMatrix[0][2] = 0;

	conf.modelParameters.transitionMatrix[1][0] = conf.modelParameters.kOff1;
	conf.modelParameters.transitionMatrix[1][1] = -conf.modelParameters.kOff1 - conf.modelParameters.kOn2;
	conf.modelParameters.transitionMatrix[1][2] = conf.modelParameters.kOn2;

	conf.modelParameters.transitionMatrix[2][0] = 0;
	conf.modelParameters.transitionMatrix[2][1] = conf.modelParameters.kOff2;
	conf.modelParameters.transitionMatrix[2][2] = -conf.modelParameters.kOff2;

	return conf;
}


//Simulation params loader
SimulationParameters load_simulationparams(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json jsonsimp = fulljson["SimulationParameters"];
	SimulationParameters simp = {};
	return assign_simulation_parameters_from_json(simp, jsonsimp);
}

//// Configuration creator new
std::vector <Configuration> load_configuration(std::string paramInputFilename, unsigned nThreads) {
	std::vector <std::string> configs = split(paramInputFilename, ";");
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
