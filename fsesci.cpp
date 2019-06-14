#include <stdexcept>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <omp.h>
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

const double E = std::exp(1.0);

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

	PotentialForce(double G_, double L_,  double sigma_)
	{
		G = G_;
		L = L_;
		powsigma = pow(sigma_, 2);
	}

	PotentialForce(double G_, double L_, double sigma_, double A_, double m_)
	{
		G = G_;
		L = L_;
		powsigma = pow(sigma_, 2);
		m = m_;
		A = A_;

		lgs = log((1.0 + m) / 2.0);
		var1 = pow(2.0, log(1.0 + m) / lgs);
		var2 = pow(2.0, log(2.0) / lgs);
	}
	
	// TODO refactor to use with period_map
	double calc(double unmodvar) const
	{
		double var = mod(unmodvar, L);
		return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);//l1d cache 4096 of doubles -> use 50% of it?
		
	}

	double asymmetric(double unmodvar) const
	{
		double x = period_map(unmodvar, L);
		double tmp1 = pow(1.0 + 2.0 * x / L, log(2.0) / lgs);

		return (log(2.0) * A * G * exp(A * (1.0 - 1.0 / (1.0 - pow((-1.0 + var1 / tmp1 ), 2)))) * tmp1 * (1.0 / tmp1 - 1.0 / var1) /
			((L + 2.0*x) * pow(-1.0 + var2 / tmp1, 2) * lgs));
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

		
	

	// rndNumbers must contain 3 * nSteps random numbers
	void advanceState(unsigned nSteps, const double* rndNumbers) {
		// configurate force object
		PotentialForce potentialForce(_mP.G, _mP.L, _mP.sigma, _mP.A, _mP.m);
		
		auto takeRandomNumber = [rndNumbers] () mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {
			
			double rnd_xMol = takeRandomNumber();
			
			double MT_Mol_force = potentialForce.asymmetric(_state.xMol);

			double next_xMol = _state.xMol + (_mP.expTime / _mP.gammaMol) * (MT_Mol_force - calculateMolspringForce(_state.xMol - _initC.xPed)) + sqrt(2*_mP.DMol*_mP.expTime) * rnd_xMol;
			
			_state.xMol   = next_xMol;
			_state.Time += _mP.expTime;
			
			_loggingBuffer.xMol   +=  _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.Time   =  _state.Time;   

		}
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_loggingBuffer);
		}
	}

	void loggingBuffertoZero() {
		_loggingBuffer.xMol = 0.0;
		_loggingBuffer.logpotentialForce = 0.0;
		_loggingBuffer.Time = 0.0;
	}

	void forcefeedbackBuffertoZero() {
		_forcefeedbackBuffer.xMol = 0.0;
		_forcefeedbackBuffer.logpotentialForce = 0.0;
		_forcefeedbackBuffer.Time = 0.0;
	
	}

private:
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
public:
	SystemState _state;
	SystemState _loggingBuffer;
	SystemState _forcefeedbackBuffer;
	const ModelParameters _mP;
	const InitialConditions _initC;
};



int main(int argc, char *argv[])
{
	if (cmdOptionExists(argv, argv + argc, "-h"))
	{
		// Do stuff
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

	int buffsize = 400'000;
	int randomsPeriter = 1;
	int stepsperbuffer = static_cast<int>(std::floor(buffsize / randomsPeriter));
	int totalsavings = 10000;//(totalsteps / iterationsbetweenSavings)
	int iterationsbetweenSavings = 15'000'000;
	int iterationsbetweenTrapsUpdate = 15'000'000;


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
		
		int macrostepMax = iterationsbetweenSavings / stepsperbuffer;

		for (int macrostep = 0; macrostep < macrostepMax; macrostep++) {
			generator1.generateNumbers();
			const auto buffData = generator1.getNumbersBuffer();

			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
				for (int taskrun = 0; taskrun < tasksperthread; taskrun++) {
					tasks[omp_get_thread_num()+ taskrun*nThreads]->advanceState(stepsperbuffer, buffData);
					
				}
			} // end of openmp section
		}
			for (const auto& task : tasks) {
				task->_forcefeedbackBuffer.xMol   += task->_loggingBuffer.xMol;
				task->_forcefeedbackBuffer.logpotentialForce += task->_loggingBuffer.logpotentialForce;

				task->_loggingBuffer.xMol   = task->_loggingBuffer.xMol / static_cast<double>(iterationsbetweenSavings);
				task->_loggingBuffer.logpotentialForce= task->_loggingBuffer.logpotentialForce / static_cast<double>(iterationsbetweenSavings);
				
				task->writeStateTolog();
				task->loggingBuffertoZero();
				
			}

		if (savedSampleIter % 10 == 0) {
			double procent = round(100 * 100 * savedSampleIter / (totalsavings)) / 100;
			std::cout << procent << "%" << std::endl;
		}
	}
}
