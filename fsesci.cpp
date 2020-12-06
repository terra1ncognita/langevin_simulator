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
		_loggedField{ loggedField }, kind{0}
	{
		auto filename = loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary";
		// print_info(filename);

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

	BinaryFileLogger(LoggerParameters loggerParams, double(MoleculeState::* loggedField), std::string coordinateName) :
		_loggedMoleculeField{ loggedField }, kind{ 1 }
	{
		auto filename = loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary";
		// print_info(filename);
		mol_num = (coordinateName.back() - '0');

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
		if (kind == 0) {
			_buffer.push_back(systemState->*_loggedField);
			
		}
		else if (kind == 1) {
			if (mol_num == 1) {
				_buffer.push_back((systemState->firstMol).*_loggedMoleculeField);
			}
			else if (mol_num == 2) {
				_buffer.push_back((systemState->secondMol).*_loggedMoleculeField);
			}
			else {
				throw std::runtime_error{ "Unknown mol_num " + std::to_string(mol_num) };
			}
		}
		else {
			throw std::runtime_error{ "Save SystemState with logger of type " + std::to_string(kind) };
		}

		if (_buffer.size() >= _buffsize) {
			flush();
		}
	}

	void save(const MoleculeState* systemState) {
		if (kind != 1) {
			throw std::runtime_error{ "Save SystemState with logger of type " + std::to_string(kind) };
		}
		_buffer.push_back(systemState->*_loggedMoleculeField);
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

	void print_info(std::string filename) {
		cout << "BinaryLogger(filename=" << filename << ", kind=" << kind << ")" << endl;
	}

	static constexpr std::size_t _buffsize = 4096 * 64 / sizeof(double);//4096 default, was 1024
	std::ofstream _file;

	double(SystemState::* _loggedField);
	double(MoleculeState::* _loggedMoleculeField);
	std::vector <double> _buffer;

	int kind = -1, mol_num;
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
	MoleculeState* state;
	const double powsigma;
	const double lgs = log(0.5), E = exp(1);
	const double var1, var2;
	const double pos = 2.0;

	PotentialForce(const ModelParameters& mp_, MoleculeState& state_) :
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
		_loggingBuffer(configuration.initialConditions.initialState),
		_forcefeedbackBuffer(configuration.initialConditions.initialState),
		expGen(1.0)
	{
		loggingBuffertoZero();
		forcefeedbackBuffertoZero();

		const auto& loggerParameters = configuration.loggerParameters;
		SystemState::iterateFields([this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});

		_state.firstMol.iterateFields([this, &loggerParameters](double(MoleculeState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});
		_state.secondMol.iterateFields([this, &loggerParameters](double(MoleculeState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		});

		_state.firstMol.expRands.assign(_mP.numStates, 0.0);
		_state.secondMol.expRands.assign(_mP.numStates, 0.0);

		_state.firstMol.livingTimes.assign(_mP.numStates, 0.0);
		_state.secondMol.livingTimes.assign(_mP.numStates, 0.0);

		fillVector(_state.firstMol.expRands);
		fillVector(_state.secondMol.expRands);

		_state.firstMol.binding = 1.0;
		if (!(_mP.bindingDynamics)) {
			_state.secondMol.binding = 1.0;
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
		PotentialForce pf1(_mP, _state.firstMol), pf2(_mP, _state.secondMol);
		std::ofstream out;
		out.open(path_prefix + "potential_" + _mP.name + "_1.txt");

		for (double x = 0; x < M_PI_2; x += 0.0001) {
			out << x << " " << pf1.well_barrier_torque(x) << endl;
		}
		out.close();

		out.open(path_prefix + "potential_" + _mP.name + "_2.txt");

		for (double x = 0; x < M_PI_2; x += 0.0001) {
			out << x << " " << pf2.well_barrier_torque(x) << endl;
		}
		out.close();
	}

	void advanceState(int nSteps, const double* const rndNumbersPointer) {
		PotentialForce potentialForce1(_mP, _state.firstMol), potentialForce2(_mP, _state.secondMol);
		
		const double* rndNumbers = rndNumbersPointer;
		auto takeRandomNumber = [rndNumbers]() mutable -> double {
			return *(rndNumbers++);
		};

		for (unsigned i = 0; i < nSteps; i++) {

			if (_mP.bindingDynamics) {
				updateState(_state.firstMol);
				updateState(_state.secondMol);
			}

			double rnd_xMT = takeRandomNumber();
			double rnd_xBeadl = takeRandomNumber();
			double rnd_xBeadr = takeRandomNumber();
			double rnd_xMol1 = takeRandomNumber();
			double rnd_phi1 = takeRandomNumber();
			double rnd_xMol2 = takeRandomNumber();
			double rnd_phi2 = takeRandomNumber();

			double MT_Mol_force_1 = potentialForce1.well_barrier_force(_state.firstMol.xMol - _state.xMT - _state.firstMol.MToffset, _state.firstMol.phi);
			double pot_torque_1 = potentialForce1.well_barrier_torque(_state.firstMol.phi);

			double MT_Mol_force_2 = potentialForce2.well_barrier_force(_state.secondMol.xMol - _state.xMT - _state.secondMol.MToffset, _state.secondMol.phi);
			double pot_torque_2 = potentialForce2.well_barrier_torque(_state.secondMol.phi);

			double molSpringForce1 = calculateMolspringForce(_state.firstMol.xMol);
			double molSpringForce2 = calculateMolspringForce(_state.secondMol.xMol);

			double MT_Mol_force = MT_Mol_force_1 + MT_Mol_force_2;

			double FmtR = calculateMTspringForce(_state.xBeadr - _state.xMT - _mP.MTlength / 2.0, 'R');
			double FmtL = calculateMTspringForce(_state.xMT - _state.xBeadl - _mP.MTlength / 2.0, 'L');
			

			double next_xMT = _state.xMT + (_sim.expTime / _mP.gammaMT)*(-FmtL + FmtR - MT_Mol_force) + sqrt(2.0*_mP.DMT*_sim.expTime) * rnd_xMT;
			double next_xBeadl = _state.xBeadl + (_sim.expTime / _mP.gammaBeadL)*((-_mP.trapstiffL)*(_state.xBeadl - _state.xTrapl) + FmtL) + sqrt(2.0*_mP.DBeadL*_sim.expTime) * rnd_xBeadl;
			double next_xBeadr = _state.xBeadr + (_sim.expTime / _mP.gammaBeadR)*(-FmtR + (-_mP.trapstiffR)*(_state.xBeadr - _state.xTrapr)) + sqrt(2.0*_mP.DBeadR*_sim.expTime) * rnd_xBeadr;
			
			double next_xMol1 = _state.firstMol.xMol + (_sim.expTime / _mP.gammaMol) * (MT_Mol_force_1 - molSpringForce1) + sqrt(2.0*_mP.DMol*_sim.expTime) * rnd_xMol1;
			double next_phi1 = _state.firstMol.phi + (_sim.expTime / _mP.rotFriction) * (-_mP.rotStiffness*(_state.firstMol.phi - _mP.iniPhi) + (_state.firstMol.binding > 0.0) * (-molSpringForce1 * _mP.molLength*sin(_state.firstMol.phi) + pot_torque_1)) + sqrt(2.0*_mP.kT*_sim.expTime / _mP.rotFriction) * rnd_phi1;

			double next_xMol2 = _state.secondMol.xMol + (_sim.expTime / _mP.gammaMol) * (MT_Mol_force_2 - molSpringForce2) + sqrt(2.0*_mP.DMol*_sim.expTime) * rnd_xMol2;
			double next_phi2 = _state.secondMol.phi + (_sim.expTime / _mP.rotFriction) * (-_mP.rotStiffness*(_state.secondMol.phi - _mP.iniPhi) + (_state.secondMol.binding > 0.0) * (-molSpringForce2 * _mP.molLength*sin(_state.secondMol.phi) + pot_torque_2)) + sqrt(2.0*_mP.kT*_sim.expTime / _mP.rotFriction) * rnd_phi2;


			if (next_phi1 < 0){
				next_phi1 = -next_phi1;
			}
			else if (next_phi1 > M_PI) {
				next_phi1 = 2 * M_PI - next_phi1;
			}

			if (next_phi2 < 0) {
				next_phi2 = -next_phi2;
			}
			else if (next_phi2 > M_PI) {
				next_phi2 = 2 * M_PI - next_phi2;
			}

			_state.xMT             = next_xMT;
			_state.xBeadl          = next_xBeadl;
			_state.xBeadr          = next_xBeadr;
			_state.firstMol.xMol   = next_xMol1;
			_state.secondMol.xMol  = next_xMol2;
			_state.firstMol.phi    = next_phi1;
			_state.secondMol.phi   = next_phi2;
			_state.Time           += _sim.expTime;

			_loggingBuffer.xMT    +=  _state.xMT;   
			_loggingBuffer.xBeadl +=  _state.xBeadl;
			_loggingBuffer.xBeadr +=  _state.xBeadr;
			_loggingBuffer.xTrapl += _state.xTrapl;
			_loggingBuffer.xTrapr += _state.xTrapr;

			_loggingBuffer.firstMol.xMol               +=  _state.firstMol.xMol;  
			_loggingBuffer.firstMol.logpotentialForce  += MT_Mol_force_1;
			_loggingBuffer.firstMol.binding            += _state.firstMol.binding;
			_loggingBuffer.firstMol.phi                += _state.firstMol.phi;
			_loggingBuffer.firstMol.potTorque          += pot_torque_1;
			_loggingBuffer.firstMol.deltaG             += _state.firstMol.deltaG;

			_loggingBuffer.secondMol.xMol              += _state.secondMol.xMol;
			_loggingBuffer.secondMol.logpotentialForce += MT_Mol_force_2;
			_loggingBuffer.secondMol.binding           += _state.secondMol.binding;
			_loggingBuffer.secondMol.phi               += _state.secondMol.phi;
			_loggingBuffer.secondMol.potTorque         += pot_torque_2;
			_loggingBuffer.secondMol.deltaG            += _state.secondMol.deltaG;
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

		_loggingBuffer.firstMol.xMol = 0.0;
		_loggingBuffer.firstMol.logpotentialForce = 0.0;
		_loggingBuffer.firstMol.binding = 0.0;
		_loggingBuffer.firstMol.phi = 0.0;
		_loggingBuffer.firstMol.potTorque = 0.0;
		_loggingBuffer.firstMol.deltaG = 0.0;

		_loggingBuffer.secondMol.xMol = 0.0;
		_loggingBuffer.secondMol.logpotentialForce = 0.0;
		_loggingBuffer.secondMol.binding = 0.0;
		_loggingBuffer.secondMol.phi = 0.0;
		_loggingBuffer.secondMol.potTorque = 0.0;
		_loggingBuffer.secondMol.deltaG = 0.0;
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
	
	void updateState(MoleculeState ms) {

		//if (ms.binding == 1.0 && (abs((ms.xMol - _state.xMT - ms.MToffset) - ms.currentWell) >= _mP.L / 2.0)) {
		//	fillVector(expRands);
		//	livingTimes.assign(_mP.numStates, 0.0);
		//	ms.binding = 0.0;
		//}
		double **transitionMatrix;
		int prev_binding = int(ms.binding);
		int	j = 0;

		if (ms.label == 1) {
			transitionMatrix = _mP.transitionMatrix1;
		}
		if (ms.label == 2) {
			transitionMatrix = _mP.transitionMatrix2;
		}

		//cout << ms.label << transitionMatrix[0][0] << transitionMatrix[0][1] << transitionMatrix[1][0] << transitionMatrix[1][1] << endl;
		cout << ms.label << endl;
		cout << &ms << endl;
		cout << &ms.expRands << endl;

		for (j = 0; j < _mP.numStates; ++j) {
			ms.livingTimes[j] += transitionMatrix[prev_binding][j] * _sim.expTime;
			cout << ms.livingTimes[j] << " ";
		}
		cout << endl;

		for (j = 0; j < _mP.numStates; ++j) {
			cout << ms.expRands[j] << ' ';
			if (ms.livingTimes[j] > ms.expRands[j] && j != prev_binding) {
				break;
			}
		}
		cout << endl;

		if (j != _mP.numStates) {
			fillVector(ms.expRands);
			ms.livingTimes.assign(_mP.numStates, 0.0);
			ms.binding = j;
		}

		//if ((prev_binding == 0.0) && (ms.binding == 1.0)) {
		//	ms.currentWell = _mP.L * floor(((ms.xMol - _state.xMT - ms.MToffset) + _mP.L / 2.0) / _mP.L);
		//}

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
};


void write_results(const std::unique_ptr<Task>& task, const SimulationParameters& sim) {
	double ibs = static_cast<double>(sim.iterationsbetweenSavings);

	task->_loggingBuffer.Time = task->_loggingBuffer.Time - sim.expTime * ibs / 2.0;

	task->_forcefeedbackBuffer.xBeadl += task->_loggingBuffer.xBeadl;
	task->_forcefeedbackBuffer.xBeadr += task->_loggingBuffer.xBeadr;
	task->_forcefeedbackBuffer.xTrapl += task->_loggingBuffer.xTrapl;
	task->_forcefeedbackBuffer.xTrapr += task->_loggingBuffer.xTrapr;

	task->_loggingBuffer.xMT    = task->_loggingBuffer.xMT / ibs;
	task->_loggingBuffer.xBeadl = task->_loggingBuffer.xBeadl / ibs;
	task->_loggingBuffer.xBeadr = task->_loggingBuffer.xBeadr / ibs;
	task->_loggingBuffer.xTrapl = task->_loggingBuffer.xTrapl / ibs;
	task->_loggingBuffer.xTrapr = task->_loggingBuffer.xTrapr / ibs;

	task->_loggingBuffer.firstMol.xMol              = task->_loggingBuffer.firstMol.xMol / ibs;
	task->_loggingBuffer.firstMol.logpotentialForce = task->_loggingBuffer.firstMol.logpotentialForce / ibs;
	task->_loggingBuffer.firstMol.binding           = task->_loggingBuffer.firstMol.binding / ibs;
	task->_loggingBuffer.firstMol.phi               = task->_loggingBuffer.firstMol.phi / ibs;
	task->_loggingBuffer.firstMol.potTorque         = task->_loggingBuffer.firstMol.potTorque / ibs;
	task->_loggingBuffer.firstMol.deltaG            = task->_loggingBuffer.firstMol.deltaG / ibs;

	task->_loggingBuffer.secondMol.xMol = task->_loggingBuffer.secondMol.xMol / ibs;
	task->_loggingBuffer.secondMol.logpotentialForce = task->_loggingBuffer.secondMol.logpotentialForce / ibs;
	task->_loggingBuffer.secondMol.binding = task->_loggingBuffer.secondMol.binding / ibs;
	task->_loggingBuffer.secondMol.phi = task->_loggingBuffer.secondMol.phi / ibs;
	task->_loggingBuffer.secondMol.potTorque = task->_loggingBuffer.secondMol.potTorque / ibs;
	task->_loggingBuffer.secondMol.deltaG = task->_loggingBuffer.secondMol.deltaG / ibs;

	task->writeStateTolog();
	task->loggingBuffertoZero();
}

void force_clamp_update(const std::unique_ptr<Task>& task, const SimulationParameters& sim) {
	double ibs = static_cast<double>(sim.iterationsbetweenSavings);

	task->_forcefeedbackBuffer.xBeadl = task->_forcefeedbackBuffer.xBeadl / ibs;
	task->_forcefeedbackBuffer.xBeadr = task->_forcefeedbackBuffer.xBeadr / ibs;
	task->_forcefeedbackBuffer.xTrapl = task->_forcefeedbackBuffer.xTrapl / ibs;
	task->_forcefeedbackBuffer.xTrapr = task->_forcefeedbackBuffer.xTrapr / ibs;

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

			// task->loggingBuffertoZero();
			// task->forcefeedbackBuffertoZero();    As long as we have buffer zeroing in Task constructor

			tasks.push_back(std::move(task));
		}
		cout << "Created list of tasks" << endl;

		cout << "Start computations..." << endl;
		for (unsigned int macrostep = 0; macrostep < sim.macrostepMax; macrostep++) {
			generator1.generateNumbers();

			#pragma omp parallel num_threads(nThreads) shared(generator1, tasks)
			{
				const double* const buffData = generator1.getNumbersBuffer();
				const auto curr_thread = omp_get_thread_num();

				for (unsigned int savedSampleIter = 0; savedSampleIter < sim.savingsPerMacrostep; savedSampleIter++) {
					std::size_t offset = savedSampleIter * sim.iterationsbetweenSavings;
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
				unsigned short minutes = static_cast<unsigned short>(floor(elapsed_seconds / 60.0));
				double seconds = elapsed_seconds - 60.0 * minutes;

				cout << "Batch " << procent << "%, total " << total_percentage << "%, elapsed " << minutes << " min " << seconds << " s" << endl;
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
