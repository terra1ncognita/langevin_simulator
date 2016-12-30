#include <stdexcept>
#include <omp.h>
//#include <stdafx.h>
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
	std::ofstream file;
	std::size_t const buffsize =4096/sizeof(double);
	std::vector <double> buffer;
	BinaryFileLogger(LoggerParameters loggerParams, std::string coordinateName) {
		
		buffer.reserve(buffsize);

		file = std::ofstream(loggerParams.filepath + loggerParams.name + "_results_"+ coordinateName+".binary" , std::ios::binary);
		if (!file) {
			throw std::runtime_error{ "the file was not created" };
		}
	}
	~BinaryFileLogger() {
		if (buffer.size()>0) {
			writeToFile();
		};
	}
	void writeToFile(){
		file.write((const char*)buffer.data(), buffer.size() * sizeof(double));
		if (!file.good()) {
			throw std::runtime_error{ "not all data was written to file" };
		};
		buffer.clear();
	}
	void save(double coordinate) {
		buffer.push_back(coordinate);
		if (buffer.size() == buffsize) {
			writeToFile();
		}
	}

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
		return (G*(-L / 2.0 + var)) / (pow(E, pow(L - 2.0 * var, 2) / (8.0*powsigma))*powsigma);
		
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
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xMT"));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xBeadl"));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xBeadr"));
			loggersvector.push_back(std::make_unique<BinaryFileLogger>(confs.at(it).loggerParameters, "xMol"));
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
#pragma omp parallel num_threads(nThreads) shared(buffData, confs, loggersvector)
				{

					int threadid = omp_get_thread_num();

					const LoggerParameters lP  = confs.at(threadid).loggerParameters;
					const ModelParameters mP   = confs.at(threadid).modelParameters;
					const InitialConditions initC = confs.at(threadid).initialConditions;
					DynamicCoordinates dynC = confs.at(threadid).dynamicCoordinates;

					// configurate force object
					PotentialForce potentialForce;
					potentialForce.E = E;
					potentialForce.G = mP.G;
					potentialForce.L = mP.L;
					potentialForce.powsigma = pow(mP.sigma, 2.0);
					
					for (int iter = 0; iter < 900'000/3; iter ++) {

						const double MT_Mol_force = potentialForce.calc(dynC.xMol - dynC.xMT);

						const double next_xMT = dynC.xMT + (sP.expTime / mP.gammaMT)*(((-mP.MTstiffL)*(dynC.xMT - dynC.xBeadl)) + (mP.MTstiffR*(dynC.xBeadr - dynC.xMT)) - (MT_Mol_force));
						const double next_xBeadl = dynC.xBeadl + (sP.expTime / mP.gammaBead)*(((-mP.trapstiff)*(dynC.xBeadl - initC.xTrapl)) + (mP.MTstiffL*(dynC.xMT - dynC.xBeadl))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter]);
						const double next_xBeadr = dynC.xBeadr + (sP.expTime / mP.gammaBead)*(((-mP.MTstiffR)*(dynC.xBeadr - dynC.xMT)) + ((-mP.trapstiff)*(dynC.xBeadr - initC.xTrapr))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter + (900'000/3)]);
						const double next_xMol = dynC.xMol + (sP.expTime / mP.gammaMol) *(MT_Mol_force + mP.molstiff*(initC.xPed - dynC.xMol)) + sqrt(2.0*mP.DMol*sP.expTime) *(buffData[iter + (2*900'000 / 3)]);

						dynC.xMT = next_xMT;
						dynC.xBeadl = next_xBeadl;
						dynC.xBeadr = next_xBeadr;
						dynC.xMol = next_xMol;
					}
					confs.at(threadid).dynamicCoordinates = dynC;
					for (int it = 0; it < confs.size(); it++) {
						loggersvector.at(4*it)->save(dynC.xMT);
						loggersvector.at(4*it+1)->save(dynC.xBeadl);
						loggersvector.at(4*it+2)->save(dynC.xBeadr);
						loggersvector.at(4*it+3)->save(dynC.xMol);
					}

					
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