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


///
double mod(const double& a, const double& N)
{
	return a - N*floor(a / N); //return in range [0, N)
}

class PotentialForce 
{
public:
	double G=0.0;
	double L=0.0;
	double E=0.0; 
	double powsigma=0.0;
	
	double calc(const double& unmodvar) const
	{
		double var = mod(unmodvar, L);
		return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);//l1d cache 4096 of doubles -> use 50% of it?		
	}
};

class Task
{
public:
	Task(const Configuration& configuration):
		
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
		
	double calculateMolspringForce(const double& extensionInput) {
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

	double calculateMTspringForce(const double& extensionInput, const char& side) {
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


	void advanceState(const unsigned& nSteps, const double* rndNumbers) {
		// configurate force object
		PotentialForce potentialForce;
		potentialForce.E = E;
		potentialForce.G = _mP.G;
		potentialForce.L = _mP.L;
		potentialForce.powsigma = pow(_mP.sigma, 2.0);

		auto takeRandomNumber = [rndNumbers] () mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {
			
			//double rnd_xMT = takeRandomNumber();
			double rnd_xMol = takeRandomNumber();
			
			double MT_Mol_force = potentialForce.calc(_state.xMol/* - _state.xMT*/);

			//double next_xMT = _state.xMT + (_mP.expTime / _mP.gammaMT)*(-calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2, 'L') + calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength/2, 'R') - (MT_Mol_force))+ sqrt(2*_mP.DMT*_mP.expTime) * rnd_xMT;
			double next_xMol = _state.xMol + (_mP.expTime / _mP.gammaMol) *(MT_Mol_force /*- calculateMolspringForce(_state.xMol - _initC.xPed)*/) + sqrt(2*_mP.DMol*_mP.expTime) * rnd_xMol;
			
			//_state.xMT    = next_xMT;
			_state.xMol   = next_xMol;
			_state.Time += _mP.expTime;


			//_loggingBuffer.xMT     +=  _state.xMT;
			_loggingBuffer.xMol    +=  _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.Time    =  _state.Time;   

		}
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_loggingBuffer);
			//logger->save(&_forcefeedbackBuffer);
		}
	}

	void loggingBuffertoZero() {
		//_loggingBuffer.xMT = 0.0;
		_loggingBuffer.xMol = 0.0;
		_loggingBuffer.logpotentialForce = 0.0;
		_loggingBuffer.Time = 0.0;
	}

	void forcefeedbackBuffertoZero() {
		//_forcefeedbackBuffer.xMT = 0.0;
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
	
	const auto configurations = load_configuration(inputparamfile,nThreads);

	std::vector<std::unique_ptr<Task>> tasks;
	for (const auto& configuration : configurations) {
		auto task = std::make_unique<Task>( configuration);
		tasks.push_back(std::move(task));
	}

	
	const int buffsize = 500'000;
	const int randomsPeriter = 1;
	const int stepsperbuffer = static_cast<int>(std::floor(buffsize / randomsPeriter));
	const int totalsavings = 7000;//(totalsteps / iterationsbetweenSavings)
	const int iterationsbetweenSavings = 15'000'000;
	const int iterationsbetweenTrapsUpdate = 15'000'000;

	const unsigned trapsUpdateTest = iterationsbetweenTrapsUpdate / iterationsbetweenSavings;
	const int macrostepMax = iterationsbetweenSavings / stepsperbuffer;

	const int tasksperthread = tasks.size() / nThreads;
	std::cout << tasksperthread << " task(s) per thread" << std::endl;


	if (iterationsbetweenSavings % stepsperbuffer != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenSavings/stepsperbuffer is integer" };
	}
	if (iterationsbetweenTrapsUpdate %  iterationsbetweenSavings != 0) {
		throw std::runtime_error{ "Please check that iterationsbetweenTrapsUpdate/iterationsbetweenSavings is integer" };
	}
	
	std::vector<MklGaussianParallelGenerator> generators(nThreads, { 0.0, 1.0, buffsize, 1});
	std::vector<const double*>buffData(nThreads);

	for (const auto& task : tasks) {
		task->loggingBuffertoZero();
		task->forcefeedbackBuffertoZero();
	}

	for (int savedSampleIter = 0; savedSampleIter < totalsavings; savedSampleIter++) {
		
		for (int macrostep = 0; macrostep < macrostepMax; macrostep++) {
			
			for (int i = 0; i < nThreads; ++i) {
				generators[i].generateNumbers();
				buffData[i] = generators[i].getNumbersBuffer();
			}

			#pragma omp parallel num_threads(nThreads) shared(buffData,tasks)
			{
				// Split tasks into consequtive sets for each thread, e.g. 0:[0 1], 1:[2 3]

				for (int taskrun = 0; taskrun < tasksperthread; taskrun++) {
				tasks[omp_get_thread_num() * tasksperthread + taskrun]->advanceState(stepsperbuffer, buffData[omp_get_thread_num()]);
				}
			} // end of openmp section

		}

		for (const auto& task : tasks) {
			//task->_forcefeedbackBuffer.xMT    += task->_loggingBuffer.xMT;
			task->_forcefeedbackBuffer.xMol   += task->_loggingBuffer.xMol;
			task->_forcefeedbackBuffer.logpotentialForce += task->_loggingBuffer.logpotentialForce;
			//task->_forcefeedbackBuffer.Time   += task->_loggingBuffer.Time;

			//task->_loggingBuffer.xMT    = task->_loggingBuffer.xMT / static_cast<double>(iterationsbetweenSavings);
			task->_loggingBuffer.xMol   = task->_loggingBuffer.xMol / static_cast<double>(iterationsbetweenSavings);
			task->_loggingBuffer.logpotentialForce= task->_loggingBuffer.logpotentialForce / static_cast<double>(iterationsbetweenSavings);
			//task->_loggingBuffer.Time   = task->_loggingBuffer.Time / static_cast<double>(iterationsbetweenSavings);

			task->writeStateTolog();
			task->loggingBuffertoZero();
				
		}
		

		if (savedSampleIter % 10 == 0) {
			double procent = round(100 * 100 * savedSampleIter / (totalsavings)) / 100;
			std::cout << procent << "%" << std::endl;
		}
	}
}