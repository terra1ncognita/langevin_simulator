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

static constexpr unsigned nThreads = 5;

// Initialize global constants
std::string inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\Langevien2_New\\Release\\config_debug.json";
const double E = std::exp(1.0);
const double kBoltz= 1.38064852e-5;// (*pN um *)


/// ToDo try to use ofstream rawwrite

class BinaryFileLogger 
{
public:
	BinaryFileLogger(LoggerParameters loggerParams, std::string coordinateName, double (SystemState::* loggedField)):
		_file{ loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary", std::ios::binary },
		_loggedField{ loggedField }
	{
		if (!_file) {
			throw std::runtime_error{ "the file was not created" };
		}
		_buffer.reserve(_buffsize);
	}
	~BinaryFileLogger() {
		writeToFile();
	}
	void writeToFile(){
		if (_buffer.empty()) {
			return;
		}
		_file.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(double));
		if (!_file.good()) {
			throw std::runtime_error{ "not all data was written to file" };
		};
		_buffer.clear();
	}
	void save(const SystemState* systemState) {
		_buffer.push_back(systemState->*_loggedField);
		if (_buffer.size() == _buffsize) {
			writeToFile();
		}
	}

private:
	static constexpr std::size_t _buffsize = 4096 / sizeof(double);
	std::ofstream _file;
	double(SystemState::* _loggedField);
	std::vector <double> _buffer;
};


///
double mod(double a, double N)
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
	
	double calc(double unmodvar) const
	{
		double var = mod(unmodvar, L);
		return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);//l1d cache 4096 of doubles -> use 50% of it?
		
	}
};





int main(int argc, char *argv[])
{
	int mode = 1;//1-> regular 2-> to check timestep
	if (mode == 1) {
		if (cmdOptionExists(argv, argv + argc, "-h"))
		{
			// Do stuff
			std::cout << "Sorry users, no help donations today." << std::endl;
		}
		char * param_input_filename = getCmdOption(argv, argv + argc, "-paramsfile");
		char * output_filename = getCmdOption(argv, argv + argc, "-resultfile");

		// Create and load simulation parameters and configuration, values are taken from json file
		SimulationParameters sP = load_simulationparams(inputparamfile);
		std::vector <Configuration> confs = load_configuration(inputparamfile);
		// Check if number of configurations correspond to predefined threads number
		if (confs.size() != nThreads) { throw std::runtime_error{ "Please check the number of configurations in json corresponds to the number of threads implemented in simulator" }; }
		
		// Create loggers for each configuration and define dynamic coordinates to be logged
		std::vector <std::unique_ptr<BinaryFileLogger>> loggersvector;
		for (int it = 0; it < confs.size(); it++) {
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xMT", &SystemState::xMT));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xBeadl", &SystemState::xBeadl));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xBeadr", &SystemState::xBeadr));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xMol", &SystemState::xMol));
		}
		///
		
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
		MklGaussianParallelGenerator generator1(0.0, 1.0, 900'000, 5);

		for (int savedstep = 0; savedstep < (100'000); savedstep++) {

			for (int macrostep = 0; macrostep < (900'000 / 900'000); macrostep++) {
				generator1.generateNumbers();
				const auto buffData = generator1.getNumbersBuffer();
#pragma omp parallel num_threads(nThreads) shared(buffData, confs)
				{

					int threadid = omp_get_thread_num();

					const ModelParameters mP   = confs.at(threadid).modelParameters;
					const InitialConditions initC = confs.at(threadid).initialConditions;
					auto currentState = confs.at(threadid).currentState;

					// configurate force object
					PotentialForce potentialForce;
					potentialForce.E = E;
					potentialForce.G = mP.G;
					potentialForce.L = mP.L;
					potentialForce.powsigma = pow(mP.sigma, 2.0);
					
					for (int iter = 0; iter < 900'000/3; iter ++) {

						const double MT_Mol_force = potentialForce.calc(currentState.xMol - currentState.xMT);

						const double next_xMT = currentState.xMT + (sP.expTime / mP.gammaMT)*(((-mP.MTstiffL)*(currentState.xMT - currentState.xBeadl)) + (mP.MTstiffR*(currentState.xBeadr - currentState.xMT)) - (MT_Mol_force));
						const double next_xBeadl = currentState.xBeadl + (sP.expTime / mP.gammaBead)*(((-mP.trapstiff)*(currentState.xBeadl - initC.xTrapl)) + (mP.MTstiffL*(currentState.xMT - currentState.xBeadl))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter]);
						const double next_xBeadr = currentState.xBeadr + (sP.expTime / mP.gammaBead)*(((-mP.MTstiffR)*(currentState.xBeadr - currentState.xMT)) + ((-mP.trapstiff)*(currentState.xBeadr - initC.xTrapr))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter + (900'000/3)]);
						const double next_xMol = currentState.xMol + (sP.expTime / mP.gammaMol) *(MT_Mol_force + mP.molstiff*(initC.xPed - currentState.xMol)) + sqrt(2.0*mP.DMol*sP.expTime) *(buffData[iter + (2*900'000 / 3)]);

						currentState.xMT = next_xMT;
						currentState.xBeadl = next_xBeadl;
						currentState.xBeadr = next_xBeadr;
						currentState.xMol = next_xMol;
					}
					confs.at(threadid).currentState = currentState;
					

					
				}//end of openmp section
				for (int it = 0; it < confs.size(); it++) {
					const auto currState = &confs.at(it).currentState;
					loggersvector.at(4 * it)->save(currState);
					loggersvector.at(4 * it + 1)->save(currState);
					loggersvector.at(4 * it + 2)->save(currState);
					loggersvector.at(4 * it + 3)->save(currState);
				}
			}
			
			if (savedstep % 100 == 0) {
				double procent = round(100*100 * savedstep / (100'000))/100;
				std::cout << procent << "%" << std::endl;
				std::cout << __rdtsc() << std::endl;
				//std::cout << nst << std::endl;
			}
		}
				
		//////////////////////
		////////////////////
		//////////////////////
	}
	
return 0;
}

//std::chrono
//_rdtscp  in #immintrin.h

//		(*/),(+-)