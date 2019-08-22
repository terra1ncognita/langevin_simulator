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
#include <math.h>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include <memory>
#include <immintrin.h>

#include <mkl_vsl.h>

#include "json.hpp"

#include "library.h"
#include "configuration.h"
#include "configuration_loader.h"
#include "mkl_gaussian_parallel_generator.h"

#include <fstream>

//static constexpr unsigned nThreads = 2;

// Initialize global constants
//std::string inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\Langevien2_New\\x64\\Release\\config_debug_trap50_noT.json";

/// ToDo try to use ofstream rawwrite

class BinaryFileLogger 
{
public:
	BinaryFileLogger(LoggerParameters loggerParams, double (SystemState::* loggedField), std::string coordinateName):
		_file{ loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary", std::ios::binary },
		_loggedField{ loggedField }
	{
		if (!_file) {
			throw std::runtime_error{ "the file was not created" };
		}
		_buffer.reserve(_buffsize);
	}
	~BinaryFileLogger() {
		flush();
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

	static constexpr std::size_t _buffsize = 4096 / sizeof(double);//4096 default, was 1024
	std::ofstream _file;
	double(SystemState::* _loggedField);
	std::vector <double> _buffer;
};


/// Old finction for symmetric well
double mod(double a, double N)
{
	return a - N*floor(a / N); //return in range [0, N)
}

// New mapping R->[-L/2, L/2], 0->0 for asymmetric well
double period_map(double y, double L)
{
	return y - floor(y / L + 0.5)* L;
}

class PotentialForce 
{
public:
	double G = 0.0, G2 = 0.0;
	double L = 0.0;
	double E = exp(1.0); 
	double powsigma = 0.0;
	double m = 0.0;
	double A = 0.0;
	double lgs = log(0.5);
	double var1 = pow(2.0, log(1.0 + m) / lgs);
	double var2 = pow(2.0, log(2.0) / lgs);
	double binding = 0.0;
	SystemState* state;

	PotentialForce(double G_, double G2_, double L_,  double sigma_, SystemState& state_)
	{
		G = G_;
		G2 = G2_;
		L = L_;
		powsigma = pow(sigma_, 2);
		state = &state_;
	}

	PotentialForce(double G_, double G2_, double L_, double sigma_, double A_, double m_, SystemState& state)
	{
		G = G_;
		G2 = G2_;
		L = L_;
		powsigma = pow(sigma_, 2);
		m = m_;
		A = A_;
		binding = state.binding;

		lgs = log((1.0 + m) / 2.0);
		var1 = pow(2.0, log(1.0 + m) / lgs);
		var2 = pow(2.0, log(2.0) / lgs);
	}
	
	double calc(double unmodvar) const
	{
		double var = period_map(unmodvar, L);
		if (state->binding == 0.0) {
			return 0.0;
		}
		else if (state->binding == 1.0) {
			return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);//l1d cache 4096 of doubles -> use 50% of it?
		}
		else if (state->binding == 2.0) {
			return (G2*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);
		}
	}

	double asymmetric(double unmodvar) const
	{
		if (binding == 0.0) {
			return 0.0;
		}
		else if (binding == 1.0) {
			double x = period_map(unmodvar, L);
			double tmp1 = pow(1.0 + 2.0 * x / L, log(2.0) / lgs);

			return (log(2.0) * A * G * exp(A * (1.0 - 1.0 / (1.0 - pow((-1.0 + var1 / tmp1), 2)))) * tmp1 * (1.0 / tmp1 - 1.0 / var1) /
				((L + 2.0*x) * pow(-1.0 + var2 / tmp1, 2) * lgs));
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
		const auto& loggerParameters = configuration.loggerParameters;
		SystemState::iterateFields([this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});
		fillVector(expRands);
	}

	double calculateMolspringForce(double extensionInput) {
		double extension = fabs(extensionInput);
		int sign = std::signbit(extensionInput);

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
		int sign = std::signbit(extensionInput);
		//(extensionInput > 0) - (extensionInput < 0);

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

	// rndNumbers must contain 3 * nSteps random numbers
	void advanceState(int nSteps, const double* rndNumbers) {
		// configurate force object
		//PotentialForce potentialForce(_mP.G, _mP.L, _mP.sigma, _mP.A, _mP.m, _state.binding);
		PotentialForce potentialForce(_mP.G, _mP.G2, _mP.L, _mP.sigma, _state);
		
		auto takeRandomNumber = [rndNumbers] () mutable -> double {
			return *(rndNumbers++);
		};

		for (int i = 0; i < nSteps; i++) {
			
			double rnd_xMT = takeRandomNumber();
			double rnd_xBeadl = takeRandomNumber();
			double rnd_xBeadr = takeRandomNumber();
			double rnd_xMol = takeRandomNumber();
			
			double MT_Mol_force = potentialForce.calc(_state.xMol - _state.xMT);
			//double MT_Mol_force = potentialForce.asymmetric(_state.xMol - _state.xMT);

			double next_xMT = _state.xMT + (_sim.expTime / _mP.gammaMT)*(-calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2.0, 'L') + calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength/2.0, 'R') - MT_Mol_force)+ sqrt(2.0*_mP.DMT*_sim.expTime) * rnd_xMT;
			double next_xBeadl = _state.xBeadl + (_sim.expTime / _mP.gammaBeadL)*(((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl)) + calculateMTspringForce(_state.xMT - _state.xBeadl -  _mP.MTlength / 2.0, 'L')) + sqrt(2.0*_mP.DBeadL*_sim.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_sim.expTime / _mP.gammaBeadR)*(-calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2.0, 'R') + ((-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr))) + sqrt(2.0*_mP.DBeadR*_sim.expTime) * rnd_xBeadr;
			double next_xMol = _state.xMol + (_sim.expTime / _mP.gammaMol) * (MT_Mol_force - calculateMolspringForce(_state.xMol - _initC.xPed)) + sqrt(2.0*_mP.DMol*_sim.expTime) * rnd_xMol;

			_state.xMT    = next_xMT;
			_state.xBeadl = next_xBeadl;
			_state.xBeadr = next_xBeadr;
			_state.xMol   = next_xMol;
			_state.Time += _sim.expTime;

			updateState();

			_loggingBuffer.xMT    +=  _state.xMT;   
			_loggingBuffer.xBeadl +=  _state.xBeadl;
			_loggingBuffer.xBeadr +=  _state.xBeadr;
			_loggingBuffer.xTrapl += _state.xTrapl;
			_loggingBuffer.xTrapr += _state.xTrapr;
			_loggingBuffer.xMol   +=  _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.Time   =  _state.Time;   
			_loggingBuffer.binding = _state.binding;
		}
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
	}

	void forcefeedbackBuffertoZero() {
		_forcefeedbackBuffer.xMT = 0.0;
		_forcefeedbackBuffer.xBeadl = 0.0;
		_forcefeedbackBuffer.xBeadr = 0.0;
		_forcefeedbackBuffer.xTrapl = 0.0;
		_forcefeedbackBuffer.xTrapr = 0.0;
		_forcefeedbackBuffer.xMol = 0.0;
		_forcefeedbackBuffer.logpotentialForce = 0.0;
		_forcefeedbackBuffer.Time = 0.0;
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
			if (livingTimes[j] > expRands[j] && j != prev_binding) {
				fillVector(expRands);
				livingTimes.assign(_mP.numStates, 0.0);
				_state.binding = j;
				break;
			}
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
	SystemState _loggingBuffer;
	SystemState _forcefeedbackBuffer;
	const SimulationParameters _sim;
	const ModelParameters _mP;
	const InitialConditions _initC;
	ExponentialGenerator expGen;
	std::vector<double> expRands;
	std::vector<double> livingTimes;
};


int main(int argc, char *argv[])
{
	if (cmdOptionExists(argv, argv + argc, "-h"))
	{
		std::cout << "Sorry users, no help donations today." << std::endl;
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

	std::vector<std::unique_ptr<Task>> tasks;
	for (const auto& configuration : configurations) {
		auto task = std::make_unique<Task>(configuration);
		tasks.push_back(std::move(task));
	}

	MklGaussianParallelGenerator generator1(0.0, 1.0, sim.buffsize, 4);

	int tasksperthread = tasks.size() / nThreads;
	std::cout << tasksperthread << std::endl;

	for (const auto& task : tasks) {
		task->loggingBuffertoZero();
		task->forcefeedbackBuffertoZero();
	}

	for (int savedSampleIter = 0; savedSampleIter < sim.totalsavings; savedSampleIter++) {

		/*for (const auto& task : tasks) {
			task->updateBinding(sim.iterationsbetweenSavings);
		}*/

		for (int macrostep = 0; macrostep < sim.macrostepMax; macrostep++) {
			generator1.generateNumbers();
			const auto buffData = generator1.getNumbersBuffer();

			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
				for (int taskrun = 0; taskrun < tasksperthread; taskrun++) {
					tasks[omp_get_thread_num()+ taskrun*nThreads]->advanceState(sim.stepsperbuffer, buffData);
				}
			} // end of openmp section
		}

		for (const auto& task : tasks) {
			task->_forcefeedbackBuffer.xMT    += task->_loggingBuffer.xMT;
			task->_forcefeedbackBuffer.xBeadl += task->_loggingBuffer.xBeadl;
			task->_forcefeedbackBuffer.xBeadr += task->_loggingBuffer.xBeadr;
			task->_forcefeedbackBuffer.xTrapl += task->_loggingBuffer.xTrapl;
			task->_forcefeedbackBuffer.xTrapr += task->_loggingBuffer.xTrapr;
			task->_forcefeedbackBuffer.xMol   += task->_loggingBuffer.xMol;
			task->_forcefeedbackBuffer.logpotentialForce += task->_loggingBuffer.logpotentialForce;
			//task->_forcefeedbackBuffer.Time   += task->_loggingBuffer.Time;

			task->_loggingBuffer.xMT    = task->_loggingBuffer.xMT / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.xBeadl = task->_loggingBuffer.xBeadl / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.xBeadr = task->_loggingBuffer.xBeadr / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.xTrapl = task->_loggingBuffer.xTrapl / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.xTrapr = task->_loggingBuffer.xTrapr / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.xMol   = task->_loggingBuffer.xMol / static_cast<double>(sim.iterationsbetweenSavings);
			task->_loggingBuffer.logpotentialForce= task->_loggingBuffer.logpotentialForce / static_cast<double>(sim.iterationsbetweenSavings);
			//task->_loggingBuffer.Time   = task->_loggingBuffer.Time / static_cast<double>(iterationsbetweenSavings);

			task->writeStateTolog();
			task->loggingBuffertoZero();
		}
		
		if ((savedSampleIter % sim.trapsUpdateTest) == 0) {
			for (const auto& task : tasks) {

				task->_forcefeedbackBuffer.xMT    = task->_forcefeedbackBuffer.xMT / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.xBeadl = task->_forcefeedbackBuffer.xBeadl / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.xBeadr = task->_forcefeedbackBuffer.xBeadr / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapl / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapr / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.xMol   = task->_forcefeedbackBuffer.xMol / static_cast<double>(sim.iterationsbetweenTrapsUpdate);
				task->_forcefeedbackBuffer.logpotentialForce = task->_forcefeedbackBuffer.logpotentialForce / static_cast<double>(sim.iterationsbetweenTrapsUpdate);

				//task->_forcefeedbackBuffer.Time   = task->_forcefeedbackBuffer.Time / static_cast<double>(iterationsbetweenTrapsUpdate);
				//task->writeStateTolog();

				int tmpDirection = task->_forcefeedbackBuffer.direction;

				if (tmpDirection == 1.0)
				{
					//moving to the right, leading bead right, trailing bead left, positive X increment
					if (task->_forcefeedbackBuffer.xBeadr >= task->_initC.initialState.xBeadr +  task->_mP.DmblMoveAmplitude) {
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
					if (task->_forcefeedbackBuffer.xBeadl <= task->_initC.initialState.xBeadl -task->_mP.DmblMoveAmplitude) {
						task->_forcefeedbackBuffer.direction = 1.0;
						task->_forcefeedbackBuffer.xTrapr = (task->_initC.initialState.xTrapr - task->_initC.initialState.xBeadr) + task->_forcefeedbackBuffer.xBeadr + (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffR);
						task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapr - task->_initC.initialState.xTrapr + task->_initC.initialState.xTrapl;
					}
					else {
						task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xBeadl + (task->_initC.initialState.xTrapl - task->_initC.initialState.xBeadl)  - (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffL);
						task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapl + (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
					}
				}
					
				task->_state.xTrapl = task->_forcefeedbackBuffer.xTrapl;
				task->_state.xTrapr = task->_forcefeedbackBuffer.xTrapr;

				task->_loggingBuffer.direction = task->_forcefeedbackBuffer.direction;
				task->forcefeedbackBuffertoZero();
			}
		}

		if (savedSampleIter % 10 == 0) {
			double procent = round(100 * 100 * savedSampleIter / sim.totalsavings) / 100;
			std::cout << procent << "%" << std::endl;
		}
	}
}
