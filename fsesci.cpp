#include <stdexcept>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string>
#include <iostream>
# include <cstdlib>
# include <iomanip>
# include <omp.h>
#include <cmath>
#include <vector>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>
#include <map>
#include <random>
#include <memory>
#include <immintrin.h>

#include <mkl_vsl.h>

#include "json.hpp"

#include "library.h"
#include "configuration.h"
#include "configuration_loader.h"
#include "mkl_gaussian_parallel_generator.h"

#include <fstream>
#include <chrono>

using std::cout;
using std::cerr;
using std::endl;

// TODO try to use ofstream rawwrite

class BinaryFileLogger 
{
public:
	BinaryFileLogger(LoggerParameters loggerParams, double (SystemState::* loggedField), std::string coordinateName):
		_loggedField{ loggedField }
	{
		auto filename = loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary";

		if (_file.is_open()) {
			_file.close();
		}

		try {
			_file.open(filename, std::ios::binary);
		}
		catch (const std::exception & ex) {
			cerr << ex.what() << endl;
			throw;
		}
		if (_file.bad()) {
			std::string err_msg = "the file " + filename + " was not created";
			cerr << err_msg << endl;
			throw std::runtime_error{ err_msg };
		}
		if (!_file.is_open()) {
			std::string err_msg = "the file " + filename + "  was not opened";
			cerr << err_msg << endl;
			throw std::runtime_error{ err_msg };
		}

		_buffer.reserve(_buffsize);
	}

	~BinaryFileLogger() {
		flush();
		_file.close();
	}

	void save(const SystemState* systemState) {
		_buffer.push_back(systemState->*_loggedField);
		if (_buffer.size() >= _buffsize) {
			flush();
		}
	}

private:
	void flush() {
		if (_buffer.empty()) {
			return;
		}
		_file.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(double));
		if (!_file.good()) {
			throw std::runtime_error{ "not all data was written to file" };
		};
		_buffer.clear();
	}

	static constexpr std::size_t _buffsize = 4096 * 64 / sizeof(double);//4096 default, was 1024
	std::ofstream _file;
	double(SystemState::* _loggedField);
	std::vector <double> _buffer;
};

// New mapping R->[-L/2, L/2], 0->0 for asymmetric well
double period_map(double y, double L)
{
	return y - floor(y / L + 0.5) * L;
}


class PotentialForce
{
public:
	const ModelParameters* mp;
	SystemState* state;
	const double powsigma;
	const double lgs = log(0.5), E = exp(1);
	const double var1, var2;
	const double pos = 2.0;


	PotentialForce(const ModelParameters& mp_, SystemState& state_) :
		powsigma ( pow(mp_.sigma, 2) ),
		var1 ( pow(2.0, log(1.0 + mp_.m) / lgs) ),
		var2 ( pow(2.0, log(2.0) / lgs) )
	{
		mp = &mp_;
		state = &state_;
	}

	
	double calc(double unmodvar) const
	{
		double var = period_map(unmodvar, mp->L);
		if (state->binding == 0.0) {
			return 0.0;
		}
		else if (state->binding == 1.0) {
			return (mp->G * var / powsigma) * pow(E, -pow(var, 2) / (2.0*powsigma));//l1d cache 4096 of doubles -> use 50% of it?
		}
		/*
		else if (state->binding == 2.0) {
			return (mp->G2 * var / powsigma) * pow(E, -pow(var, 2) / (2.0*powsigma));
		}
		*/
	}
	
	double asymmetric(double unmodvar) const
	{
		if (state->binding == 0.0) {
			return 0.0;
		}
		else if (state->binding == 1.0) {
			double x = period_map(unmodvar, mp->L);
			double tmp1 = pow(1.0 + 2.0 * x / mp->L, log(2.0) / lgs);

			return (log(2.0) * mp->A * mp->G * exp(mp->A * (1.0 - 1.0 / (1.0 - pow((-1.0 + var1 / tmp1), 2)))) * tmp1 * (1.0 / tmp1 - 1.0 / var1) /
				((mp->L + 2.0*x) * pow(-1.0 + var2 / tmp1, 2) * lgs));
		}
	}

	double calc_potential_torque(double angle) const
	{
		if (angle >= M_PI_2) {
			return 0.0;
		}
		return -(pow(mp->domainsDistance, 2) * exp(2 - 2 * mp->domainsDistance * sin(angle) / mp->rotWellWidth) *
			2 * mp->rotWellDepth * (mp->rotWellWidth - mp->domainsDistance * sin(angle)) * sin(2 * angle)) / pow(mp->domainsDistance, 3);
	}

	double two_domains(double unmodvar, double angle) const
	{
		double var = period_map(unmodvar, mp->L);
		double deltaG = 0.0;

		if (state->binding == 0.0) {
			return 0.0;
		}

		if (angle < M_PI_2) {
			deltaG = 2 * mp->rotWellDepth * pow((mp->domainsDistance * sin(angle) / mp->rotWellWidth), 2) * exp(2 * (1 - mp->domainsDistance * sin(angle) / mp->rotWellWidth));
		}

		state->deltaG = deltaG;

		if (state->binding == 1.0) {
			return ((mp->G + deltaG) * var / powsigma) * pow(E, -pow(var, 2) / (2.0*powsigma));//l1d cache 4096 of doubles -> use 50% of it?
		}
	}

	double morze(double x, double r0, double depth) const
	{
		return depth * pow((x / r0), 2) * exp(2 * (1 - x / r0));
	}

	double morze_angle_derivative(double angle, double r0, double d, double depth) const
		/*
		Assume relationship x =  mp->domainsDistance * sin(angle)
		*/
	{
		return -depth * (pow(d, 2) * exp(2 - 2 * d * sin(angle) / r0) *
			(r0 - d * sin(angle)) * sin(2 * angle)) / pow(r0, 3);
	}

	double well_barrier_torque(double angle) const
	{
		if (angle >= M_PI_2) {
			return 0.0;
		}
		return morze_angle_derivative(angle, mp->rotWellWidth, mp->domainsDistance, mp->rotWellDepth) +
			morze_angle_derivative(angle, pos * mp->rotWellWidth, mp->domainsDistance, -mp->B * mp->rotWellDepth);
	}

	double well_barrier_force(double unmodvar, double angle) const
	{
		double var = period_map(unmodvar, mp->L);
		double deltaG = 0.0;

		if (state->binding == 0.0) {
			return 0.0;
		}

		if (angle < M_PI_2) {
			deltaG = morze(mp->domainsDistance * sin(angle), mp->rotWellWidth, mp->rotWellDepth) +
				morze(mp->domainsDistance * sin(angle), pos * mp->rotWellWidth, -mp->B * mp->rotWellDepth);
		}

		state->deltaG = deltaG;

		if (state->binding == 1.0) {
			return ((mp->G + deltaG) * var / powsigma) * pow(E, -pow(var, 2) / (2.0*powsigma));//l1d cache 4096 of doubles -> use 50% of it?
		}
	}
};

class ExponentialGenerator {
public:
	ExponentialGenerator(double lambda) :
		expDist(lambda),
		re(std::random_device{}())
	{
	}

	double operator()() {
		return expDist(re);
	}

private:
	std::default_random_engine re;
	std::exponential_distribution<double> expDist;
};


class Task
{
public:
	Task(const Configuration& configuration) :
		_sim(configuration.simulationParameters),
		_mP(configuration.modelParameters),
		_initC(configuration.initialConditions),
		_state(configuration.initialConditions.initialState),
		expGen(1.0),
		expRands(_mP.numStates),
		livingTimes(_mP.numStates, 0.0)
	{
		loggingBuffertoZero();
		forcefeedbackBuffertoZero();

		const auto& loggerParameters = configuration.loggerParameters;
		SystemState::iterateFields([this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});

		fillVector(expRands);

		if (!(_mP.bindingDynamics)) {
			_state.binding = 1.0;
		}
	}

	double calculateMolspringForce(double extensionInput) {
		double extension = fabs(extensionInput);
		int sign = (extensionInput > 0) - (extensionInput < 0);

		if (extension <= _mP.molStiffBoundary) {
			return sign * _mP.molStiffWeakSlope*extension;
		}
		else
		{
			return sign * (_mP.molStiffStrongSlope*extension + _mP.molStiffBoundary*(_mP.molStiffWeakSlope - _mP.molStiffStrongSlope));
		}

	}

	double calculateMTspringForce(double extensionInput, char side) {
		double extension = fabs(extensionInput)*2.0;
		int sign = (extensionInput > 0) - (extensionInput < 0);

		if (side == 'L')
		{
			if (extension <= _mP.MTstiffWeakBoundaryL)
			{
				return _mP.MTstiffWeakSlopeL*extension*sign;
			}
			else
			{
				if (extension > _mP.MTstiffWeakBoundaryL && extension < _mP.MTstiffStrongBoundaryL)
				{
					return (_mP.MTstiffParabolicAL*pow(extension, 2) + _mP.MTstiffParabolicBL*extension + _mP.MTstiffParabolicCL)*sign;
				}
				else
				{
					return (_mP.MTstiffStrongSlopeL*extension + _mP.MTstiffStrongIntersectL)*sign;
				}
			}
		}
		if (side == 'R')
		{
			if (extension <= _mP.MTstiffWeakBoundaryR)
			{
				return _mP.MTstiffWeakSlopeR*extension*sign;
			}
			else
			{
				if (extension > _mP.MTstiffWeakBoundaryR && extension < _mP.MTstiffStrongBoundaryR)
				{
					return (_mP.MTstiffParabolicAR*pow(extension, 2) + _mP.MTstiffParabolicBR*extension + _mP.MTstiffParabolicCR)*sign;
				}
				else
				{
					return (_mP.MTstiffStrongSlopeR*extension + _mP.MTstiffStrongIntersectR)*sign;
				}
			}
		}

	}

	void log_well_torque(std::string path_prefix) {
		PotentialForce pf(_mP, _state);
		std::ofstream out;
		out.open(path_prefix + "potential_" + _mP.name + ".txt");

		for (double x = 0; x < M_PI_2; x += 0.0001) {
			out << x << " " << pf.well_barrier_torque(x) << endl;
		}
		out.close();
	}

	void advanceState(int nSteps, const double* rndNumbers) {
		PotentialForce potentialForce(_mP, _state);
		
		auto takeRandomNumber = [rndNumbers]() mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {

			if (_mP.bindingDynamics) {
				updateState();
			}

			double rnd_xMT = takeRandomNumber();
			double rnd_xBeadl = takeRandomNumber();
			double rnd_xBeadr = takeRandomNumber();
			double rnd_xMol = takeRandomNumber();
			double rnd_phi = takeRandomNumber();

			//double MT_Mol_force = potentialForce.calc(_state.xMol - _state.xMT);
			//double MT_Mol_force = potentialForce.asymmetric(_state.xMol - _state.xMT);
			double MT_Mol_force = potentialForce.well_barrier_force(_state.xMol - _state.xMT, _state.phi);
			double pot_torque = potentialForce.well_barrier_torque(_state.phi);

			double FmtR = calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2.0, 'R');
			double FmtL = calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2.0, 'L');
			double molSpringForce = calculateMolspringForce(_state.xMol);

			double next_xMT = _state.xMT + (_sim.expTime / _mP.gammaMT)*(-FmtL + FmtR - MT_Mol_force) + sqrt(2.0*_mP.DMT*_sim.expTime) * rnd_xMT;
			double next_xBeadl = _state.xBeadl + (_sim.expTime / _mP.gammaBeadL)*((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl) + FmtL) + sqrt(2.0*_mP.DBeadL*_sim.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_sim.expTime / _mP.gammaBeadR)*(-FmtR + (-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr)) + sqrt(2.0*_mP.DBeadR*_sim.expTime) * rnd_xBeadr;
			double next_xMol = _state.xMol + (_sim.expTime / _mP.gammaMol) * (MT_Mol_force - molSpringForce) + sqrt(2.0*_mP.DMol*_sim.expTime) * rnd_xMol;
			double next_phi = _state.phi + (_sim.expTime / _mP.rotFriction) * (-_mP.rotStiffness*(_state.phi - _mP.iniPhi) + (_state.binding > 0.0) * (-molSpringForce * _mP.molLength*sin(_state.phi) + pot_torque)) + sqrt(2.0*_mP.kT*_sim.expTime / _mP.rotFriction) * rnd_phi;

			if (next_phi < 0){
				next_phi = -next_phi;
			}
			else if (next_phi > M_PI) {
				next_phi = 2 * M_PI - next_phi;
			}

			_state.xMT    = next_xMT;
			_state.xBeadl = next_xBeadl;
			_state.xBeadr = next_xBeadr;
			_state.xMol   = next_xMol;
			_state.phi    = next_phi;
			_state.Time += _sim.expTime;

			_loggingBuffer.xMT    +=  _state.xMT;   
			_loggingBuffer.xBeadl +=  _state.xBeadl;
			_loggingBuffer.xBeadr +=  _state.xBeadr;
			_loggingBuffer.xTrapl += _state.xTrapl;
			_loggingBuffer.xTrapr += _state.xTrapr;
			_loggingBuffer.xMol   +=  _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.binding += _state.binding;
			_loggingBuffer.phi += _state.phi;
			_loggingBuffer.potTorque += pot_torque;
			_loggingBuffer.deltaG += _state.deltaG;
		}
		_loggingBuffer.Time = _state.Time;
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_loggingBuffer);
		}
	}

	void loggingBuffertoZero() {
		_loggingBuffer.xMT = 0.0;
		_loggingBuffer.xBeadl = 0.0;
		_loggingBuffer.xBeadr = 0.0;
		_loggingBuffer.xTrapl = 0.0;
		_loggingBuffer.xTrapr = 0.0;
		_loggingBuffer.xMol = 0.0;
		_loggingBuffer.logpotentialForce = 0.0;
		_loggingBuffer.Time = 0.0;
		_loggingBuffer.binding = 0.0;
		_loggingBuffer.phi = 0.0;
		_loggingBuffer.potTorque = 0.0;
		_loggingBuffer.deltaG = 0.0;
	}

	void forcefeedbackBuffertoZero() {
		_forcefeedbackBuffer.xBeadl = 0.0;
		_forcefeedbackBuffer.xBeadr = 0.0;
		_forcefeedbackBuffer.xTrapl = 0.0;
		_forcefeedbackBuffer.xTrapr = 0.0;
	}
	
	void fillVector(std::vector<double>& rnds) {
		std::generate(begin(rnds), end(rnds), expGen);
	}
	
	void updateState() {

		if (_state.binding == 1.0 && (abs((_state.xMol - _state.xMT) - _state.currentWell) >= _mP.L / 2.0)) {
			fillVector(expRands);
			livingTimes.assign(_mP.numStates, 0.0);
			_state.binding = 0.0;
		}

		int prev_binding = int(_state.binding);
		int	j = 0;

		for (j = 0; j < _mP.numStates; ++j) {
			livingTimes[j] += _mP.transitionMatrix[prev_binding][j] * _sim.expTime;
		}

		for (j = 0; j < _mP.numStates; ++j) {
			if (livingTimes[j] > expRands[j] && j != prev_binding) {
				break;
			}
		}

		if (j != _mP.numStates) {
			fillVector(expRands);
			livingTimes.assign(_mP.numStates, 0.0);
			_state.binding = j;
		}

		if ((prev_binding == 0.0) && (_state.binding == 1.0)) {
			_state.currentWell = _mP.L * floor(((_state.xMol - _state.xMT) + _mP.L / 2.0) / _mP.L);
		}

		return;
	}

private:
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
public:
	SystemState _state;
	SystemState _loggingBuffer, _forcefeedbackBuffer;
	const SimulationParameters _sim;
	const ModelParameters _mP;
	const InitialConditions _initC;
	ExponentialGenerator expGen;
	std::vector<double> expRands;
	std::vector<double> livingTimes;
};


void write_results(const std::unique_ptr<Task>& task, const SimulationParameters& sim) {
	task->_forcefeedbackBuffer.xBeadl += task->_loggingBuffer.xBeadl;
	task->_forcefeedbackBuffer.xBeadr += task->_loggingBuffer.xBeadr;
	task->_forcefeedbackBuffer.xTrapl += task->_loggingBuffer.xTrapl;
	task->_forcefeedbackBuffer.xTrapr += task->_loggingBuffer.xTrapr;

	task->_loggingBuffer.xMT = task->_loggingBuffer.xMT / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.xBeadl = task->_loggingBuffer.xBeadl / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.xBeadr = task->_loggingBuffer.xBeadr / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.xTrapl = task->_loggingBuffer.xTrapl / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.xTrapr = task->_loggingBuffer.xTrapr / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.xMol = task->_loggingBuffer.xMol / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.logpotentialForce = task->_loggingBuffer.logpotentialForce / static_cast<double>(sim.iterationsbetweenSavings);

	task->_loggingBuffer.Time = task->_loggingBuffer.Time - sim.expTime * static_cast<double>(sim.iterationsbetweenSavings) / 2;

	task->_loggingBuffer.binding = task->_loggingBuffer.binding / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.phi = task->_loggingBuffer.phi / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.potTorque = task->_loggingBuffer.potTorque / static_cast<double>(sim.iterationsbetweenSavings);
	task->_loggingBuffer.deltaG = task->_loggingBuffer.deltaG / static_cast<double>(sim.iterationsbetweenSavings);

	task->writeStateTolog();
	task->loggingBuffertoZero();
}

void force_clamp_update(const std::unique_ptr<Task>& task, const SimulationParameters& sim) {
	task->_forcefeedbackBuffer.xBeadl = task->_forcefeedbackBuffer.xBeadl / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
	task->_forcefeedbackBuffer.xBeadr = task->_forcefeedbackBuffer.xBeadr / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
	task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapl / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
	task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapr / static_cast<double>(sim.iterationsbetweenTrapsUpdate);

	int tmpDirection = task->_forcefeedbackBuffer.direction;

	if (tmpDirection == 1.0)
	{
		//moving to the right, leading bead right, trailing bead left, positive X increment
		if (task->_forcefeedbackBuffer.xBeadr >= task->_initC.initialState.xBeadr + task->_mP.DmblMoveAmplitude) {
			task->_forcefeedbackBuffer.direction = -1.0;
			task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xBeadl + (task->_initC.initialState.xTrapl - task->_initC.initialState.xBeadl) - (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffL);
			task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapl + (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
		}
		else {
			task->_forcefeedbackBuffer.xTrapr = (task->_initC.initialState.xTrapr - task->_initC.initialState.xBeadr) + task->_forcefeedbackBuffer.xBeadr + (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffR);
			task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapr - (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
		}
	}
	if (tmpDirection == -1.0)
	{
		//moving to the left, leading bead left, trailing bead right, negative X increment
		if (task->_forcefeedbackBuffer.xBeadl <= task->_initC.initialState.xBeadl - task->_mP.DmblMoveAmplitude) {
			task->_forcefeedbackBuffer.direction = 1.0;
			task->_forcefeedbackBuffer.xTrapr = (task->_initC.initialState.xTrapr - task->_initC.initialState.xBeadr) + task->_forcefeedbackBuffer.xBeadr + (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffR);
			task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapr - task->_initC.initialState.xTrapr + task->_initC.initialState.xTrapl;
		}
		else {
			task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xBeadl + (task->_initC.initialState.xTrapl - task->_initC.initialState.xBeadl) - (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffL);
			task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapl + (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
		}
	}

	task->_state.xTrapl = task->_forcefeedbackBuffer.xTrapl;
	task->_state.xTrapr = task->_forcefeedbackBuffer.xTrapr;

	task->_loggingBuffer.direction = task->_forcefeedbackBuffer.direction;
	task->forcefeedbackBuffertoZero();
}


int main(int argc, char *argv[])
{
	auto start_main = std::chrono::system_clock::now();

	if (cmdOptionExists(argv, argv + argc, "-h"))
	{
		cout << "Sorry users, no help donations today." << endl;
	}
	char* param_input_filename = getCmdOption(argv, argv + argc, "-paramsfile");
	char* output_filename = getCmdOption(argv, argv + argc, "-resultfile");
	char* charnThreads = getCmdOption(argv, argv + argc, "-nthreads");

	std::string inputparamfile(param_input_filename);
	unsigned nThreads = std::stoi(charnThreads);
	
	// Create and load simulation parameters and configuration, values are taken from json file
	const auto configurations = load_configuration(inputparamfile, nThreads);

	SimulationParameters sim;
	if (cmdOptionExists(argv, argv + argc, "-taskfile")) {
		char* sim_filename = getCmdOption(argv, argv + argc, "-taskfile");
		sim = load_simulationparams(std::string(sim_filename));
	}
	else {
		sim = SimulationParameters();
	}
	cout << "Loaded simultaion parameters" << endl;

	MklGaussianParallelGenerator generator1(0.0, 1.0, sim.buffsize, 4);
	cout << "Created random numbers generator" << endl;

	int tasksperthread = configurations.size() / nThreads;
	cout << "Total number of batches is " << tasksperthread << endl;
	cout << std::setprecision(2) << std::fixed << endl;

	for (int batch_id = 0; batch_id < tasksperthread; batch_id++) {
		// Loop of task batches. One batch fully laods all threads with exactly one task.

		cout << "Start batch #" << batch_id << endl;
		auto start_batch = std::chrono::system_clock::now();

		std::vector<Configuration> batch;
		batch.reserve(nThreads);

		for (int conf_id = batch_id * nThreads; conf_id < (batch_id + 1) * nThreads; conf_id++) {
			batch.push_back(std::move(configurations[conf_id]));
		}

		cout << "Contains " << batch.size() << " simulations: ";
		for (const auto& conf:batch) {
			cout << conf.loggerParameters.name << " ";
		}
		cout << endl;
		
		std::vector<std::unique_ptr<Task>> tasks;
		tasks.reserve(batch.size());

		for (const auto& configuration : batch) {
			std::unique_ptr<Task> task = std::make_unique<Task>(configuration);

			task->loggingBuffertoZero();
			task->forcefeedbackBuffertoZero();

			tasks.push_back(std::move(task));
		}

		for (int macrostep = 0; macrostep < sim.macrostepMax; macrostep++) {
			generator1.generateNumbers();
			const auto buffData = generator1.getNumbersBuffer();

			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
				for (int savedSampleIter = 0; savedSampleIter < sim.savingsPerMacrostep; savedSampleIter++) {
					tasks[omp_get_thread_num()]->advanceState(sim.iterationsbetweenSavings, buffData);
					write_results(tasks[omp_get_thread_num()], sim);

					if ((savedSampleIter % sim.trapsUpdateTest) == 0) {
						force_clamp_update(tasks[omp_get_thread_num()], sim);
					}
				}
			} // end of openmp section

			int counter = macrostep + 1;
			if (counter % 20 == 0) {
				double procent = 100.0 * static_cast<double>(counter) / sim.macrostepMax;
				double total_percentage = (100.0 * batch_id + procent) / tasksperthread;

				auto curr = std::chrono::system_clock::now();
				std::chrono::duration<double> dt = (curr - start_batch);
				double elapsed_seconds = dt.count();
				int minutes = static_cast<int>(floor(elapsed_seconds / 60.0));
				double seconds = elapsed_seconds - 60.0 * minutes;

				cout << "Batch " << procent << "%, total " << total_percentage << "%, elapsed " << minutes << " min " << seconds << " s" << endl;
			}
		}

		auto curr = std::chrono::system_clock::now();
		std::chrono::duration<double> dt = (curr - start_batch);
		double elapsed_seconds = dt.count();
		int minutes = static_cast<int>(floor(elapsed_seconds / 60));
		double seconds = elapsed_seconds - 60 * minutes;

		cout << endl << "Finished batch #" << batch_id << " in " << minutes << " min " << seconds << " s" << endl;

		dt = (curr - start_main);
		elapsed_seconds = dt.count();
		minutes = static_cast<int>(floor(elapsed_seconds / 60));
		seconds = elapsed_seconds - 60 * minutes;
		cout << "Elapsed " << minutes << " min " << seconds << " s from program start" << endl << endl;
	}

	auto curr = std::chrono::system_clock::now();
	std::chrono::duration<double> dt = (curr - start_main);
	double elapsed_seconds = dt.count();
	int minutes = static_cast<int>(floor(elapsed_seconds / 60));
	double seconds = elapsed_seconds - 60 * minutes;
	cout << endl << "All simulations finished in " << minutes << " min " << seconds << " s from program start" << endl;

	return 0;
}
