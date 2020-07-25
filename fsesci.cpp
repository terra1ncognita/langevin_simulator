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

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>
#include <map>
#include <random>
#include <memory>
#include <immintrin.h>

#include <mkl_vsl.h>

//#include "json.hpp"

#include "library.h"
#include "configuration.h"
#include "configuration_loader.h"
#include "mkl_gaussian_parallel_generator.h"
#include "binding_state.h"

#include <fstream>
#include <chrono>

//#include "binding_state.h"

using std::cout;
using std::cerr;
using std::endl;

// TODO try to use ofstream rawwrite

template<typename T>
class BinaryFileLogger 
{
public:
	BinaryFileLogger() {};
	BinaryFileLogger(LoggerParameters loggerParams, T (SystemState::* loggedField), std::string coordinateName):
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
		_file.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(T));
		if (!_file.good()) {
			throw std::runtime_error{ "not all data was written to file" };
		};
		_buffer.clear();
	}

	static constexpr std::size_t _buffsize = 4096 * 64 / sizeof(T);//4096 default, was 1024
	std::ofstream _file;
	T(SystemState::* _loggedField);
	std::vector <T> _buffer;
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

	double x_rot_min = -1, x_rot_max = -1;

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
		if (state->molecular_state == BindingState::Free) {
			return 0.0;
		}
		else if (state->molecular_state == BindingState::FirstSiteBound) {
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
		if (state->molecular_state == BindingState::Free) {
			return 0.0;
		}
		else if (state->molecular_state & BindingState::FirstSiteBound) {
			double x = period_map(unmodvar, mp->L);
			double tmp1 = pow(1.0 + 2.0 * x / mp->L, log(2.0) / lgs);

			return (log(2.0) * mp->A * mp->G * exp(mp->A * (1.0 - 1.0 / (1.0 - pow((-1.0 + var1 / tmp1), 2)))) * tmp1 * (1.0 / tmp1 - 1.0 / var1) /
				((mp->L + 2.0*x) * pow(-1.0 + var2 / tmp1, 2) * lgs));
		}
	}

	double morze(double x, double r0, double depth) const
	{
		return depth * pow((x / r0), 2) * exp(2 * (1 - x / r0));
	}

	double morze_derivative(double x, double r0, double depth) const
	{
		/*
		Returns morze'(x)
		*/
		return 2.0 * exp(2.0 - 2.0 * x / r0) * depth * (r0 - x) * x / pow(r0, 3);
	}

	double well_barrier_derivative(double x) const
	{
		return morze_derivative(x, mp->rotWellWidth, mp->rotWellDepth) +
			morze_derivative(x, pos * mp->rotWellWidth, -mp->B * mp->rotWellDepth);
	}

	double morze_angle_neg_derivative(double angle, double r0, double d, double depth) const
		/*
		Returns -morze'(angle) assuming relationship x = d * sin(angle).
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
		return morze_angle_neg_derivative(angle, mp->rotWellWidth, mp->domainsDistance, mp->rotWellDepth) +
			morze_angle_neg_derivative(angle, pos * mp->rotWellWidth, mp->domainsDistance, -mp->B * mp->rotWellDepth);
	}

	double well_barrier_force(double unmodvar, double angle) const
	{
		double var = period_map(unmodvar, mp->L);
		double deltaG = 0.0;

		if (state->molecular_state == BindingState::Free) {
			return 0.0;
		}

		if (angle < M_PI_2) {
			deltaG = morze(mp->domainsDistance * sin(angle), mp->rotWellWidth, mp->rotWellDepth) +
				morze(mp->domainsDistance * sin(angle), pos * mp->rotWellWidth, -mp->B * mp->rotWellDepth);
		}

		state->deltaG = deltaG;

		if (state->molecular_state & BindingState::FirstSiteBound) {
			return ((mp->G + deltaG) * var / powsigma) * pow(E, -pow(var, 2) / (2.0*powsigma)); //l1d cache 4096 of doubles -> use 50% of it?
		}
	}

	double solve_rot_well_newton(double g, double x0, double tol) const 
	{
		int max_iter = 1000;
		double x=x0, x_next = x0;

		for (int i=0; i < max_iter; ++i)
		{
			x = x_next;
			double f_val = morze(x, mp->rotWellWidth, mp->rotWellDepth) + morze(x, pos * mp->rotWellWidth, -mp->B * mp->rotWellDepth) - g;
			if (-tol < f_val && f_val < tol)
			{
				break;
			}

			x_next = x - (f_val) 
				/ (morze_derivative(x, mp->rotWellWidth, mp->rotWellDepth) + morze_derivative(x, pos * mp->rotWellWidth, -mp->B * mp->rotWellDepth));
		}
		return x_next;
	}

	double secant_solver_well_derivative(double x0, double x1, double tol) const
	{
		int max_iter = 1000;
		double x, x_next, f_prev, f_val, tmp;

		x = x0;
		x_next = x1;

		for (int i = 0; i < max_iter; ++i)
		{
			f_prev = well_barrier_derivative(x);
			f_val = well_barrier_derivative(x_next);

			if (-tol < f_val && f_val < tol)
			{
				break;
			}

			tmp = x_next;
			x_next = x_next - f_val * (x_next - x) / (f_val - f_prev);
			x = tmp;
		}
		return x_next;
	}

	double get_rot_well_min_max(double rmax0, double rmin0, double eps, double tol)
	{
		/* Perturb initial values for r with eps to satisfy secant method 
		Approach max from the left, min from the right = start from right well wall
		*/
		x_rot_max = secant_solver_well_derivative(rmax0, rmax0 + eps, tol);
		x_rot_min = secant_solver_well_derivative(rmin0 + eps, rmin0, tol);
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
		livingTimes(_mP.numStates, 0.0),
		_bindingEventLogger(configuration.loggerParameters, &SystemState::lastBinding, "bindingEvent"),
		_releaseEventLogger(configuration.loggerParameters, &SystemState::lastRelease, "releaseEvent"),
		_eventTimeLogger(configuration.loggerParameters, &SystemState::Time, "eventTime")
	{
		loggingBuffertoZero();

		const auto& loggerParameters = configuration.loggerParameters;
		SystemState::iterateFields([this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger<double>>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});

		fillVector(expRands);

		if (!(_mP.bindingDynamics)) {
			_state.binding = 1.0;
			_state.molecular_state |= BindingState::FirstSiteBound;
		}
	}

	/*double calculateMolspringForce(double extensionInput) {
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
*/

	void log_well_torque(std::string path_prefix) {
		PotentialForce pf(_mP, _state);
		std::ofstream out;
		out.open(path_prefix + "potential_" + _mP.name + ".txt");

		for (double x = 0; x < M_PI_2; x += 0.0001) {
			out << x << " " << pf.well_barrier_torque(x) << endl;
		}
		out.close();
	}

	BindingEvent check_binding_event(double pot_torque) {
		bitmask::bitmask<BindingState> change = _state.molecular_state;

		// if ((_mP.domainsDistance * sin(_state.phi) > 2 * _mP.rotWellWidth) && (pot_torque < 0)) {  // TODO check this fucking condition, I hate it
		// 	change &= ~BindingState::SecondSiteBound;
		// }
		// else {
		// 	change |= BindingState::SecondSiteBound;
		// }

		if (_mP.domainsDistance * sin(_state.phi) > rotWellReleaseBoundary) {
			change &= ~BindingState::SecondSiteBound;
		}
		if (_mP.domainsDistance * sin(_state.phi) < rotWellBindingBoundary) {
			change |= BindingState::SecondSiteBound;
		}


		// if (abs(period_map(_state.xMol - _state.xMT, _mP.L)) > 3 * _mP.sigma) {
		// 	change &= ~BindingState::FirstSiteBound;
		// }
		// else {
		// 	change |= BindingState::FirstSiteBound;
		// }

		if (abs(period_map(_state.xMol - _state.xMT, _mP.L)) < _mP.sigma / 2.0) {
			change |= BindingState::FirstSiteBound;
			if (abs(_state.xMol - _state.xMT - _state.currentWell) > _mP.L / 2.0) {
				change |= BindingState::NewWell;
			}
		}
		if (abs(period_map(_state.xMol - _state.xMT, _mP.L)) > 3.0 * _mP.sigma) {
			change &= ~BindingState::FirstSiteBound;
		}

		auto event = detect_changes(_state.molecular_state, change);
		_state.molecular_state = change;

		return event;
	}

	void init_well_boundaries() {
		PotentialForce potentialForce(_mP, _state);

		potentialForce.get_rot_well_min_max(1.99 * _mP.rotWellWidth, 1.01 * _mP.rotWellWidth, 0.01 * _mP.rotWellWidth, 1e-4);

		const double g_rel = 0.8825;
		double max_height = potentialForce.morze(potentialForce.x_rot_max, _mP.rotWellWidth, _mP.rotWellDepth) + potentialForce.morze(potentialForce.x_rot_max, 2 * _mP.rotWellWidth, -_mP.B * _mP.rotWellDepth);
		double min_height = potentialForce.morze(potentialForce.x_rot_min, _mP.rotWellWidth, _mP.rotWellDepth) + potentialForce.morze(potentialForce.x_rot_min, 2 * _mP.rotWellWidth, -_mP.B * _mP.rotWellDepth);

		double starting_point;
		if (_mP.B == 0.0) {
			max_height = 0;
			starting_point = potentialForce.x_rot_min + _mP.rotWellWidth / 2.0;
		}
		else {
			starting_point = (potentialForce.x_rot_max + potentialForce.x_rot_min) / 2.0;
		}

		rotWellBindingBoundary = potentialForce.solve_rot_well_newton((1 - g_rel) * max_height + g_rel * min_height, starting_point, 1e-4);
		rotWellReleaseBoundary = potentialForce.x_rot_max;
	}

	void advanceState(int nSteps, const double* rndNumbers) {
		PotentialForce potentialForce(_mP, _state);
		bool bound_flg = true;

		auto takeRandomNumber = [rndNumbers]() mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {

			if (_mP.bindingDynamics) {
				updateState();
				bound_flg = _state.binding > 0;
			}

			double rnd_xMol = takeRandomNumber();
			double rnd_phi = takeRandomNumber();

			double MT_Mol_force = potentialForce.well_barrier_force(_state.xMol - _state.xMT, _state.phi);
			double pot_torque = potentialForce.well_barrier_torque(_state.phi);

			BindingEvent ev = check_binding_event(pot_torque);
			if (ev.any_change()) {
				_state.lastBinding = ev.binding;
				_state.lastRelease = ev.release;
				_bindingEventLogger.save(&_state);
				_releaseEventLogger.save(&_state);
				_eventTimeLogger.save(&_state);

				if (ev.binding & BindingState::NewWell) {
					_state.currentWell = _mP.L * floor(((_state.xMol - _state.xMT) + _mP.L / 2.0) / _mP.L);
					_state.molecular_state &= ~BindingState::NewWell;
				}
			}

			double next_xMol = _state.xMol + (_sim.expTime / _mP.gammaMol) * (MT_Mol_force) + sqrt(2.0*_mP.DMol*_sim.expTime) * rnd_xMol;
			double next_phi = _state.phi + (_sim.expTime / _mP.rotFriction) * (-_mP.rotStiffness*(_state.phi - _mP.iniPhi) + bound_flg * (pot_torque)) + sqrt(2.0*_mP.kT*_sim.expTime / _mP.rotFriction) * rnd_phi;

			if (next_phi < 0){
				next_phi = -next_phi;
			}
			else if (next_phi > M_PI) {
				next_phi = 2 * M_PI - next_phi;
			}

			_state.xMol  = next_xMol;
			_state.phi   = next_phi;
			_state.Time += _sim.expTime;

			_loggingBuffer.xMol              += _state.xMol;  
			_loggingBuffer.logpotentialForce += MT_Mol_force;
			_loggingBuffer.binding           += _state.binding;
			_loggingBuffer.phi               += _state.phi;
			_loggingBuffer.potTorque         += pot_torque;
			_loggingBuffer.deltaG            += _state.deltaG;
		}
		_loggingBuffer.Time = _state.Time;
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
		_loggingBuffer.binding = 0.0;
		_loggingBuffer.phi = 0.0;
		_loggingBuffer.potTorque = 0.0;
		_loggingBuffer.deltaG = 0.0;
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
	std::vector<std::unique_ptr<BinaryFileLogger<double>>> _loggers;
	BinaryFileLogger<bitmask::bitmask<BindingState>> _bindingEventLogger, _releaseEventLogger;
	BinaryFileLogger<double> _eventTimeLogger;

	double rotWellReleaseBoundary, rotWellBindingBoundary;
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

	MklGaussianParallelGenerator generator1(0.0, 1.0, sim.buffsize * nThreads, 4);
	cout << "Created random numbers generator" << endl;

	int tasksperthread = configurations.size() / nThreads;
	cout << "Total number of batches is " << tasksperthread << endl;
	cout << std::setprecision(2) << std::fixed << endl;

	for (int batch_id = 0; batch_id < tasksperthread; batch_id++) {
		// Loop over task batches. One batch fully laods all threads with exactly one task.

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
			task->init_well_boundaries();
			tasks.push_back(std::move(task));
		}

		for (int macrostep = 0; macrostep < sim.macrostepMax; macrostep++) {
			generator1.generateNumbers();
			const auto buffData = generator1.getNumbersBuffer();

			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
				for (int savedSampleIter = 0; savedSampleIter < sim.savingsPerMacrostep; savedSampleIter++) {
					tasks[omp_get_thread_num()]->advanceState(sim.iterationsbetweenSavings, buffData + sim.buffsize * omp_get_thread_num());
					write_results(tasks[omp_get_thread_num()], sim);
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
