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
const double E = std::exp(1.0);
//const double kBoltz= 1.38064852e-5; //pN um


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

	static constexpr std::size_t _buffsize = 1024 / sizeof(double);//4096 default
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
	double x = y + L / 2.0;
	return x - floor(x / L)* L - L / 2.0;
}

class PotentialForce 
{
public:
	double G = 0.0;
	double L = 0.0;
	double E = exp(1.0); 
	double powsigma = 0.0;
	double m = 0.0;
	double A = 0.0;
	double lgs = log(0.5);
	double var1 = pow(2.0, log(1.0 + m) / lgs);
	double var2 = pow(2.0, log(2.0) / lgs);
	double binding = 0.0;

	PotentialForce(double G_, double L_,  double sigma_, double binding_)
	{
		G = G_;
		L = L_;
		powsigma = pow(sigma_, 2);
		binding = binding_;
	}

	PotentialForce(double G_, double L_, double sigma_, double A_, double m_, double binding_)
	{
		G = G_;
		L = L_;
		powsigma = pow(sigma_, 2);
		m = m_;
		A = A_;
		binding = binding_;

		lgs = log((1.0 + m) / 2.0);
		var1 = pow(2.0, log(1.0 + m) / lgs);
		var2 = pow(2.0, log(2.0) / lgs);
	}
	
	// TODO refactor to use with period_map
	double calc(double unmodvar) const
	{
		if (binding == 0.0) {
			return 0.0;
		}
		else if (binding == 1.0) {
			double var = mod(unmodvar, L);
			return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);//l1d cache 4096 of doubles -> use 50% of it?
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

class Task
{
public:
	Task( const Configuration& configuration):
		
		_mP( configuration.modelParameters ),
		_initC( configuration.initialConditions ),
		_state( configuration.initialConditions.initialState )
		
	{
		const auto& loggerParameters = configuration.loggerParameters;
		SystemState::iterateFields([this, &loggerParameters] (double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});
	}
		
	double calculateMolspringForce(double extensionInput) {
		double extension = fabs(extensionInput);
		int sign = (extensionInput > 0) - (extensionInput < 0);

		if (extension <= _mP.molStiffBoundary) {
			return sign*_mP.molStiffWeakSlope*extension;
		}
		else
		{
			return sign*(_mP.molStiffStrongSlope*extension + _mP.molStiffBoundary*(_mP.molStiffWeakSlope-_mP.molStiffStrongSlope));
		}
		 
	}

	double calculateMTspringForce(double extensionInput, char side) {
		//double extensionNm = fabs(extension * 1000);
		//return (extension/fabs(extension))*(0.0062*extensionNm+(1.529*pow(10,-6)*(pow(extensionNm,2))) + (2.72*pow(10,-7)*pow(extensionNm,3)));
		double extension=fabs(extensionInput)*2.0;
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
					return (_mP.MTstiffParabolicAL*pow(extension,2) + _mP.MTstiffParabolicBL*extension + _mP.MTstiffParabolicCL)*sign;
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
	void advanceState(unsigned nSteps, const double* rndNumbers) {
		// configurate force object
		//PotentialForce potentialForce(_mP.G, _mP.L, _mP.sigma, _mP.A, _mP.m, _state.binding);
		PotentialForce potentialForce(_mP.G, _mP.L, _mP.sigma, _state.binding);

		auto takeRandomNumber = [rndNumbers] () mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {
			
			double rnd_xMT = takeRandomNumber();
			double rnd_xBeadl = takeRandomNumber();
			double rnd_xBeadr = takeRandomNumber();
			double rnd_xMol = takeRandomNumber();
			
			//double MT_Mol_force = potentialForce.calc(_state.xMol - _state.xMT);

			double MT_Mol_force = potentialForce.asymmetric(_state.xMol - _state.xMT);

			double next_xMT = _state.xMT + (_mP.expTime / _mP.gammaMT)*(-calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2, 'L') + calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength/2, 'R') - MT_Mol_force)+ sqrt(2*_mP.DMT*_mP.expTime) * rnd_xMT;
			double next_xBeadl = _state.xBeadl + (_mP.expTime / _mP.gammaBeadL)*(((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl)) + calculateMTspringForce(_state.xMT - _state.xBeadl -  _mP.MTlength / 2, 'L')) + sqrt(2*_mP.DBeadL*_mP.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_mP.expTime / _mP.gammaBeadR)*(-calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2, 'R') + ((-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr))) + sqrt(2*_mP.DBeadR*_mP.expTime) * rnd_xBeadr;
			double next_xMol = _state.xMol + (_mP.expTime / _mP.gammaMol) * (MT_Mol_force - calculateMolspringForce(_state.xMol - _initC.xPed)) + sqrt(2*_mP.DMol*_mP.expTime) * rnd_xMol;
			
			//double next_vMol = (next_xMol - _state.xMol) / _mP.expTime;
			//double next_vMT = (next_xMT - _state.xMT) / _mP.expTime;

			/*double next_xBeadl = _state.xBeadl + (_mP.expTime / _mP.gammaBeadL)*(((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl)) + calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2, 'L')) + sqrt(2 * _mP.DBeadL*_mP.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_mP.expTime / _mP.gammaBeadR)*(-calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2, 'R') + ((-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr))) + sqrt(2 * _mP.DBeadR*_mP.expTime) * rnd_xBeadr;

			double next_xMT = (_mP.expTime / _mP.gammaMT)*(-calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2, 'L') + calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2, 'R')) + sqrt(2 * _mP.DMT*_mP.expTime) * rnd_xMT;
			double next_xMol = (_mP.expTime / _mP.gammaMol) *(- calculateMolspringForce(_state.xMol - _initC.xPed)) + sqrt(2 * _mP.DMol*_mP.expTime) * rnd_xMol;
			
			next_xMol = _state.xMol + (_mP.gammaQuasiviscous * _mP.gammaMT * next_xMT + (_mP.gammaMT + _mP.gammaQuasiviscous) * _mP.gammaMol * next_xMol) / (_mP.gammaMT * _mP.gammaQuasiviscous + _mP.gammaMol * (_mP.gammaMT + _mP.gammaQuasiviscous));
			next_xMT = _state.xMT + (_mP.gammaQuasiviscous * _mP.gammaMol * next_xMol + (_mP.gammaMol + _mP.gammaQuasiviscous) * _mP.gammaMT * next_xMT) / (_mP.gammaMT * _mP.gammaQuasiviscous + _mP.gammaMol * (_mP.gammaMT + _mP.gammaQuasiviscous));
			
			double next_vMol = (next_xMol - _state.xMol) / _mP.expTime;
			double next_vMT = (next_xMT - _state.xMT) / _mP.expTime;*/

			_state.xMT    = next_xMT;
			_state.xBeadl = next_xBeadl;
			_state.xBeadr = next_xBeadr;
			_state.xMol   = next_xMol;
			_state.Time += _mP.expTime;
			//_state.vMol = next_vMol;
			//_state.vMT = next_vMT;

			//double MT_Mol_force = _mP.gammaQuasiviscous * (_state.vMT - _state.vMol);

			_loggingBuffer.xMT    +=  _state.xMT;   
			_loggingBuffer.xBeadl +=  _state.xBeadl;
			_loggingBuffer.xBeadr +=  _state.xBeadr;
			_loggingBuffer.xTrapl += _state.xTrapl;
			_loggingBuffer.xTrapr += _state.xTrapr;
			_loggingBuffer.xMol   +=  _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.Time   =  _state.Time;   

		}
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_loggingBuffer);
			//logger->save(&_forcefeedbackBuffer);
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

	void updadeBinding(int steps) {
		double p;
		double _rnd = uniform_unit(re);

		if (_state.binding == 0.0) {
			p = 1 - exp(-_mP.kOn * _mP.expTime * steps);
			if (p > _rnd) {
				_state.binding == 1.0;
			}
		}
		else if (_state.binding == 1.0) {
			p = 1 - exp(-_mP.kOff * _mP.expTime * steps);
			if (p > _rnd) {
				_state.binding == 0.0;
			}
		}
	}

private:
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
public:
	SystemState _state;
	SystemState _loggingBuffer;
	SystemState _forcefeedbackBuffer;
	const ModelParameters _mP;
	const InitialConditions _initC;
	const std::uniform_real_distribution<double> uniform_unit;
	const std::default_random_engine re;
};



int main(int argc, char *argv[])
{
	if (cmdOptionExists(argv, argv + argc, "-h"))
	{
		std::cout << "Sorry users, no help donations today." << std::endl;
	}
	char * param_input_filename = getCmdOption(argv, argv + argc, "-paramsfile");
	char * output_filename = getCmdOption(argv, argv + argc, "-resultfile");
	char * charnThreads = getCmdOption(argv, argv + argc, "-nthreads");

	std::string inputparamfile;
	inputparamfile.append(param_input_filename);
	unsigned nThreads = std::stoi(charnThreads);
	
	// Create and load simulation parameters and configuration, values are taken from json file
	//const auto simulationParameters = load_simulationparams(inputparamfile);
	const auto configurations = load_configuration(inputparamfile,nThreads);

	std::vector<std::unique_ptr<Task>> tasks;
	for (const auto& configuration : configurations) {
		auto task = std::make_unique<Task>( configuration);
		tasks.push_back(std::move(task));
	}

	// Check if number of configurations correspond to predefined threads number
	/*if (configurations.size() != nThreads) {
		throw std::runtime_error{ "Please check the number of configurations in json corresponds to the number of threads implemented in simulator" };
	}
	*/
	/*

	///////////////////////////////
	//////////////////////////////
	///////////// Iterations
	///////////////////////////
	/*
	const int frequency = (int)round(sP.microsteps);
	const int totalsimsteps = (int)round((sP.microsteps)*(sP.nTotal));
	const int nsteps = (int)round(3*(totalsimsteps / frequency) / intbufferSize);
		*/


		/*
		int frequency = sP.microsteps;
		int totalsimsteps = sP.microsteps*(sP.nTotal);
		int nsteps = 3 * (totalsimsteps / frequency) / intbufferSize;
		int dd = 2;

		frequency 100'000
		totalsimsteps 100'000'000'000
		buffer 1'000'000
		randms per simstep 3

		saved step range is 10'000 -> nTotal
		microsteps is macrostep*(buffersize/3)
		*/

	int buffsize = 800'000;
	int randomsPeriter = 4;
	int stepsperbuffer = static_cast<int>(std::floor(buffsize / randomsPeriter));
	int totalsavings = 7000;//(totalsteps / iterationsbetweenSavings)
	int iterationsbetweenSavings = 15'000'000;
	int iterationsbetweenTrapsUpdate = 15'000'000;
	int macrostepMax = iterationsbetweenSavings / stepsperbuffer;

	if (iterationsbetweenSavings % stepsperbuffer != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenSavings/stepsperbuffer is integer" };
	}
	if (iterationsbetweenTrapsUpdate %  iterationsbetweenSavings != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenTrapsUpdate/iterationsbetweenSavings is integer" };
	}

	unsigned trapsUpdateTest = iterationsbetweenTrapsUpdate / iterationsbetweenSavings;

	MklGaussianParallelGenerator generator1(0.0, 1.0, buffsize, 4);

	int tasksperthread = tasks.size() / nThreads;
	std::cout << tasksperthread << std::endl;

	for (const auto& task : tasks) {
		task->loggingBuffertoZero();
		task->forcefeedbackBuffertoZero();
		
	}

	for (int savedSampleIter = 0; savedSampleIter < totalsavings; savedSampleIter++) {

		for (const auto& task : tasks) {
			task->updadeBinding(iterationsbetweenSavings);
		}

		for (int macrostep = 0; macrostep < macrostepMax; macrostep++) {
			generator1.generateNumbers();
			const auto buffData = generator1.getNumbersBuffer();

			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
					//const auto begin = __rdtsc();
				//tasks[omp_get_thread_num()]->advanceState(stepsperbuffer, buffData);
				for (int taskrun = 0; taskrun < tasksperthread; taskrun++) {
					tasks[omp_get_thread_num()+ taskrun*nThreads]->advanceState(stepsperbuffer, buffData);
					
				}
					//const auto end = __rdtsc();
					//const auto cyclesPerStep = static_cast<double>(end - begin) / static_cast<double>(std::floor(buffsize / randomsPeriter));
					//std::cout << "cyclesPerStep = " << cyclesPerStep << std::endl;
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

				task->_loggingBuffer.xMT    = task->_loggingBuffer.xMT / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.xBeadl = task->_loggingBuffer.xBeadl / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.xBeadr = task->_loggingBuffer.xBeadr / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.xTrapl = task->_loggingBuffer.xTrapl / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.xTrapr = task->_loggingBuffer.xTrapr / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.xMol   = task->_loggingBuffer.xMol / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.logpotentialForce= task->_loggingBuffer.logpotentialForce / static_cast<double>(iterationsbetweenSavings);
				//task->_loggingBuffer.Time   = task->_loggingBuffer.Time / static_cast<double>(iterationsbetweenSavings);

				task->writeStateTolog();
				task->loggingBuffertoZero();
				
			}
		
			 if ((savedSampleIter % trapsUpdateTest) == 0) {
				for (const auto& task : tasks) {

					task->_forcefeedbackBuffer.xMT    = task->_forcefeedbackBuffer.xMT / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.xBeadl = task->_forcefeedbackBuffer.xBeadl / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.xBeadr = task->_forcefeedbackBuffer.xBeadr / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapl / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapr / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.xMol   = task->_forcefeedbackBuffer.xMol / static_cast<double>(iterationsbetweenTrapsUpdate);
					task->_forcefeedbackBuffer.logpotentialForce = task->_forcefeedbackBuffer.logpotentialForce / static_cast<double>(iterationsbetweenTrapsUpdate);

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
			double procent = round(100 * 100 * savedSampleIter / (totalsavings)) / 100;
			std::cout << procent << "%" << std::endl;
			//std::cout << __rdtsc() << std::endl;
			//std::cout << nst << std::endl;
		}
	}
}

//std::chrono
//_rdtscp  in #immintrin.h

//		(*/),(+-)