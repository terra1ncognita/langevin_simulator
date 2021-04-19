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
#include "polynomial.h"
#include "mkl_gaussian_parallel_generator.h"

#include <fstream>
#include <chrono>

#include<utility>

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



class Task
{
public:
	Task(const Configuration& configuration) :
		_sim(configuration.simulationParameters),
		_mP(configuration.modelParameters),
		_initC(configuration.initialConditions),
		_state(configuration.initialConditions.initialState),
		_loggingBuffer(configuration.initialConditions.initialState),
		_forcefeedbackBuffer(configuration.initialConditions.initialState),
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

		if (!(_state.isFree)) {
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


	void advanceState(int nSteps, const double* const rndNumbersPointer) {
		const double* rndNumbers = rndNumbersPointer;
		auto takeRandomNumber = [rndNumbers]() mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {

			double rnd_xMT = takeRandomNumber();
			double rnd_xBeadl = takeRandomNumber();
			double rnd_xBeadr = takeRandomNumber();

			double FmtR = _mP.fec(_state.xBeadr - _state.xMT - _mP.MTlength / 2.0); // Right
			double FmtL = _mP.fec(_state.xMT - _state.xBeadl - _mP.MTlength / 2.0); // Left

			double molSpringForce = 0.0;
			if (!_state.isFree) {
				molSpringForce = calculateMolspringForce(_state.xMT - _state.xMTbinding);
			}

			double next_xMT = _state.xMT + (_sim.expTime / _mP.gammaMT)*(-FmtL + FmtR - molSpringForce) + sqrt(2.0*_mP.DMT*_sim.expTime) * rnd_xMT;
			double next_xBeadl = _state.xBeadl + (_sim.expTime / _mP.gammaBeadL)*((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl) + FmtL) + sqrt(2.0*_mP.DBeadL*_sim.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_sim.expTime / _mP.gammaBeadR)*(-FmtR + (-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr)) + sqrt(2.0*_mP.DBeadR*_sim.expTime) * rnd_xBeadr;
		
			_state.xMT    = next_xMT;
			_state.xBeadl = next_xBeadl;
			_state.xBeadr = next_xBeadr;
			_state.Time += _sim.expTime;

			_loggingBuffer.xMT    +=  _state.xMT;   
			_loggingBuffer.xBeadl +=  _state.xBeadl;
			_loggingBuffer.xBeadr +=  _state.xBeadr;
			_loggingBuffer.xTrapl += _state.xTrapl;
			_loggingBuffer.xTrapr += _state.xTrapr;
			_loggingBuffer.binding += _state.binding;
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
		_loggingBuffer.Time = 0.0;
		_loggingBuffer.binding = 0.0;
	}

	void forcefeedbackBuffertoZero() {
		_forcefeedbackBuffer.xBeadl = 0.0;
		_forcefeedbackBuffer.xBeadr = 0.0;
		_forcefeedbackBuffer.xTrapl = 0.0;
		_forcefeedbackBuffer.xTrapr = 0.0;
	}
	

private:
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
public:
	SystemState _state;
	SystemState _loggingBuffer, _forcefeedbackBuffer;
	const SimulationParameters _sim;
	const ModelParameters _mP;
	const InitialConditions _initC;
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

	task->_loggingBuffer.Time = task->_loggingBuffer.Time - sim.expTime * static_cast<double>(sim.iterationsbetweenSavings) / 2;

	task->_loggingBuffer.binding = task->_loggingBuffer.binding / static_cast<double>(sim.iterationsbetweenSavings);

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
		task->_forcefeedbackBuffer.xTrapr = (task->_initC.initialState.xTrapr - task->_initC.initialState.xBeadr) + task->_forcefeedbackBuffer.xBeadr + (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffR);
		task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapr - (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
	}
	else if (tmpDirection == -1.0)
	{
		//moving to the left, leading bead left, trailing bead right, negative X increment
		task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xBeadl + (task->_initC.initialState.xTrapl - task->_initC.initialState.xBeadl) - (0.5*task->_mP.movementTotalForce / task->_mP.trapstiffL);
		task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapl + (task->_initC.initialState.xTrapr - task->_initC.initialState.xTrapl);
	}

	task->_state.xTrapl = task->_forcefeedbackBuffer.xTrapl;
	task->_state.xTrapr = task->_forcefeedbackBuffer.xTrapr;

	task->_loggingBuffer.direction = task->_forcefeedbackBuffer.direction;
	task->forcefeedbackBuffertoZero();
}

struct MinSec {
	unsigned short minutes;
	double seconds;

	MinSec(double elapsed_seconds) {
		minutes = static_cast<unsigned short>(floor(elapsed_seconds / 60.0));
		seconds = elapsed_seconds - 60.0 * minutes;
	}
};

std::ostream& operator<< (std::ostream& out, MinSec obj) {
	out << obj.minutes << " min " << obj.seconds << " s";
	return out;
}


int main(int argc, char *argv[])
{
	std::ifstream in("D:\\test.json");
	json js;
	in >> js;
	CubicSpline cs(js);

	cout << cs.find_domain(-5.4) << endl;
	cout << cs.find_domain(-5.29) << endl;
	cout << cs.find_domain(0.169) << endl;
	cout << cs.find_domain(0.2) << endl;

	cout << "Boundaries " << cs.left_boundary << " " << cs.right_boundary << endl;

	cout << "Evaluate" << endl;

	for (double x = cs.left_boundary; x <= cs.right_boundary; x += 0.01) {
		cout << x << " " << cs(x) << endl;
	}

	return 0;

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

	MklGaussianParallelGenerator generator1(0.0, 1.0, sim.buffsize * nThreads, sim.rndThreads);
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
			tasks.push_back(std::move(task));
		}
		cout << "Created list of tasks" << endl;

		ForceExtensionCurve ff = tasks[0]->_mP.fec;

		for (double x = ff.left_pole(); x < ff.right_pole(); x += 0.001) {
			cout << x << " " << ff(x) << endl;
		}
		return 0;

		cout << "Start computations..." << endl;
		for (unsigned int macrostep = 0; macrostep < sim.macrostepMax; macrostep++) {
			generator1.generateNumbers();

			for (auto &task : tasks) {
				if (macrostep < sim.macrostepsFree) {
					task->_state.binding = 0.0;
					task->_state.instaBindingDynamics = false;
					task->_state.isFree = true;
				}
				else if (task->_state.isFree) {
					task->_state.instaBindingDynamics = task->_mP.bindingDynamics;
					task->_state.binding = 1.0;
					task->_state.isFree = false;
					task->_state.xMTbinding = task->_state.xMT;
				}
			}

			#pragma omp parallel num_threads(nThreads) shared(generator1, tasks)
			{
				const double* const buffData = generator1.getNumbersBuffer();
				const auto curr_thread = omp_get_thread_num();

				for (unsigned int savedSampleIter = 0; savedSampleIter < sim.savingsPerMacrostep; savedSampleIter++) {
					std::size_t offset = savedSampleIter * sim.iterationsbetweenSavings + sim.buffsize * curr_thread;
					const double* const rnd_pointer = buffData + offset;

					tasks[curr_thread]->advanceState(sim.iterationsbetweenSavings, rnd_pointer);
					write_results(tasks[curr_thread], sim);

					if ((savedSampleIter % sim.trapsUpdateTest) == 0) {
						force_clamp_update(tasks[curr_thread], sim);
					}
				}
			} // end of openmp section

			unsigned int counter = macrostep + 1;
			if (counter % 20 == 0) {
				double procent = 100.0 * static_cast<double>(counter) / sim.macrostepMax;
				double total_percentage = (100.0 * batch_id + procent) / tasksperthread;

				auto curr = std::chrono::system_clock::now();
				std::chrono::duration<double> dt = (curr - start_batch);

				double elapsed_seconds = dt.count();
				double eta_seconds = (100.0 / total_percentage - 1) * elapsed_seconds;
				MinSec elapsed(elapsed_seconds), eta(eta_seconds);

				cout << "Batch " << procent << "%, total " << total_percentage << "%, elapsed " << elapsed << ", total ETA " << eta << endl;
			}
		}

		auto curr = std::chrono::system_clock::now();
		std::chrono::duration<double> dt = (curr - start_batch);
		double elapsed_seconds = dt.count();
		unsigned short minutes = static_cast<unsigned short>(floor(elapsed_seconds / 60));
		double seconds = elapsed_seconds - 60 * minutes;

		cout << endl << "Finished batch #" << batch_id << " in " << minutes << " min " << seconds << " s" << endl;

		dt = (curr - start_main);
		elapsed_seconds = dt.count();
		minutes = static_cast<unsigned short>(floor(elapsed_seconds / 60));
		seconds = elapsed_seconds - 60 * minutes;
		cout << "Elapsed " << minutes << " min " << seconds << " s from program start" << endl << endl;
	}

	auto curr = std::chrono::system_clock::now();
	std::chrono::duration<double> dt = (curr - start_main);
	double elapsed_seconds = dt.count();
	unsigned short minutes = static_cast<unsigned short>(floor(elapsed_seconds / 60));
	double seconds = elapsed_seconds - 60 * minutes;
	cout << endl << "All simulations finished in " << minutes << " min " << seconds << " s from program start" << endl;

	return 0;
}
