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
		simp.buffsize = static_cast<std::size_t>(stoi(jsonobjsimp["buffsize"].get<std::string>()));
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
	if (!(jsonobj["InitialConditions"]["xMol1"].empty())) {
		conf.initialConditions.initialState.firstMol.xMol = stod(jsonobj["InitialConditions"]["xMol1"].get<std::string>());
		conf.currentState.firstMol.xMol = stod(jsonobj["InitialConditions"]["xMol1"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMol2"].empty())) {
		conf.initialConditions.initialState.secondMol.xMol = stod(jsonobj["InitialConditions"]["xMol2"].get<std::string>());
		conf.currentState.secondMol.xMol = stod(jsonobj["InitialConditions"]["xMol2"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMT"].empty())) {
		conf.initialConditions.initialState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
		conf.currentState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadl"].empty())) {
		conf.initialConditions.initialState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
		conf.currentState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadr"].empty())) {
		conf.initialConditions.initialState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
		conf.currentState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapl"].empty())) {
		conf.initialConditions.initialState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
		conf.currentState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.initialConditions.initialState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
		conf.currentState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["direction"].empty())) {
		conf.currentState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
		conf.initialConditions.initialState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["phi1"].empty())) {
		conf.currentState.firstMol.phi = stod(jsonobj["InitialConditions"]["phi1"].get<std::string>());
		conf.initialConditions.initialState.secondMol.phi = stod(jsonobj["InitialConditions"]["phi1"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["phi2"].empty())) {
		conf.currentState.secondMol.phi = stod(jsonobj["InitialConditions"]["phi2"].get<std::string>());
		conf.initialConditions.initialState.secondMol.phi = stod(jsonobj["InitialConditions"]["phi2"].get<std::string>());
	}

	//TODO: generalize to arbitrary Markov chain
	int nstates = conf.modelParameters.numStates;
	if (nstates != 2) {
		throw std::runtime_error{ "Current implementation is only for 2 states!" };
	}
	conf.modelParameters.transitionMatrix1 = ModelParameters::assign_rates(nstates, conf.modelParameters.kOn1, conf.modelParameters.kOff1);
	conf.modelParameters.transitionMatrix2 = ModelParameters::assign_rates(nstates, conf.modelParameters.kOn2, conf.modelParameters.kOff2);

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
