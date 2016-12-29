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

static constexpr unsigned nThreads = 5;

// Initialize global constants
std::string inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\Langevien2_New\\Release\\config_debug.json";
const double E = std::exp(1.0);
const double kBoltz= 1.38064852e-5;// (*pN um *)


/// ToDo try to use ofstream rawwrite

class BinaryFileLogger 
{
public:
	FILE* pFile0;
	FILE* pFile1;
	FILE* pFile2;
	FILE* pFile3;
	std::vector <FILE*> files;

	std::size_t const buffsize =128;
	std::vector <double> buffer0;
	std::vector <double> buffer1;
	std::vector <double> buffer2;
	std::vector <double> buffer3;
	std::vector <std::vector <double>> buffers;
	BinaryFileLogger(LoggerParameters loggerParams) {

		pFile0 = fopen((loggerParams.filepath + loggerParams.name + std::string{ "_results_xMT.binary" }).c_str(), "wb");
		pFile1 = fopen((loggerParams.filepath + loggerParams.name + std::string{ "_results_xBeadl.binary" }).c_str(), "wb");
		pFile2 = fopen((loggerParams.filepath + loggerParams.name + std::string{ "_results_xBeadr.binary" }).c_str(), "wb");
		pFile3 = fopen((loggerParams.filepath + loggerParams.name + std::string{ "_results_xMol.binary" }).c_str(), "wb");
		buffer0.reserve(buffsize);
		buffer1.reserve(buffsize);
		buffer2.reserve(buffsize);
		buffer3.reserve(buffsize);
		if (pFile0 == nullptr) {
			throw std::runtime_error{ "the file was not created" };
		}
	}
	void save(double xMT, double xBeadl, double xBeadr, double xMol) {
		buffer0.push_back(xMT);
		if (buffer0.size() == buffsize) {
			
			if (fwrite(buffer0.data(), sizeof(double), buffsize, pFile0) != buffsize) {
				throw std::runtime_error{ "not all data was written to file" };
			};
			
			buffer0.clear();
		}
		buffer1.push_back(xBeadl);
		if (buffer1.size() == buffsize) {

			if (fwrite(buffer1.data(), sizeof(double), buffsize, pFile1) != buffsize) {
				throw std::runtime_error{ "not all data was written to file" };
			};

			buffer1.clear();
		}
		buffer2.push_back(xBeadr);
		if (buffer2.size() == buffsize) {

			if (fwrite(buffer2.data(), sizeof(double), buffsize, pFile2) != buffsize) {
				throw std::runtime_error{ "not all data was written to file" };
			};

			buffer2.clear();
		}
		buffer3.push_back(xMol);
		if (buffer3.size() == buffsize) {

			if (fwrite(buffer3.data(), sizeof(double), buffsize, pFile3) != buffsize) {
				throw std::runtime_error{ "not all data was written to file" };
			};

			buffer3.clear();
		}
	}
	void closeLogger() {
		if (buffer0.size()>0) {
			fwrite(buffer0.data(), sizeof(double), buffer0.size(), pFile0);
		};
		buffer0.clear();
		if (buffer1.size()>0) {
			fwrite(buffer1.data(), sizeof(double), buffer1.size(), pFile1);
		};
		buffer1.clear();
		if (buffer2.size()>0) {
			fwrite(buffer2.data(), sizeof(double), buffer2.size(), pFile2);
		};
		buffer2.clear();
		if (buffer3.size()>0) {
			fwrite(buffer3.data(), sizeof(double), buffer3.size(), pFile3);
		};
		buffer3.clear();
		fclose(pFile0);
		fclose(pFile1);
		fclose(pFile2);
		fclose(pFile3);
	}
};
/*
BinaryFileLogger::BinaryFileLogger(LoggerParameters loggerParams) {
	
	pFile = fopen((loggerParams.filepath + loggerParams.name + std::string{ "_results.binary" }).c_str(), "wb");
	buffer.reserve(buffsize);
	if (pFile == nullptr) {
		throw std::runtime_error{ "the file was not created" };
	}
};
*/

BinaryFileLogger create_logger(LoggerParameters loggerParams) {
	BinaryFileLogger logger(loggerParams);
	return logger;
}

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

		//////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////
		////////////////// Random Number Generator  ////////////////

		VSLStreamStatePtr stream = nullptr;
		const auto creationResult = vslNewStream(&stream, VSL_BRNG_SFMT19937, 777);
		if (creationResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "can't create VSL stream" };
		}

		static constexpr std::size_t buffSize = 900000;// Maximum size should fit L3 chache //1200'000
		
		//std::size_t buffSize = 1200'000;
		std::vector<double> buff(buffSize);

		const auto generateNumbers = [&buff, stream]() {

			static constexpr std::size_t nPerThread = buffSize / nThreads;
			static_assert(buffSize == nPerThread * nThreads, "buffsize must be multiple of nThreads");

			const auto buffData = buff.data();
			const auto _stream = stream;
#pragma omp parallel num_threads(nThreads) shared(buffData, _stream)
			{
				const std::size_t begin = nPerThread * omp_get_thread_num();
				const auto generationResult = vdRngGaussian(
					VSL_RNG_METHOD_GAUSSIAN_ICDF, _stream,
					nPerThread, buffData + begin,
					0.0, 1.0
				);
				if (generationResult != VSL_STATUS_OK) {
					throw std::runtime_error{ "can't generate numbers" };
				}
			}

		};

		////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////












		//const ModelParameters mP;
		//InitialConditions iC;
		//const SimulationParameters sP;

		SimulationParameters sP = load_simulationparams(inputparamfile);
		std::vector <Configuration> conf = load_configuration(inputparamfile);

		//const auto logger = createLogger(LoggerType::BINARY_FILE, lP);
		std::vector <BinaryFileLogger> loggersvector;
		for (int it = 0; it <= 4; it++) {
			loggersvector.push_back(create_logger(conf.at(it).loggerParameters));
		}

		//		(*/),(+-)

		///////////////////////////////
		//////////////////////////////
		///////////// Iterations
		///////////////////////////
		/*
		const int frequency = (int)round(sP.microsteps);
		const int totalsimsteps = (int)round((sP.microsteps)*(sP.nTotal));
		const int nsteps = (int)round(3*(totalsimsteps / frequency) / intbufferSize);
			*/

		double xcoords[nThreads][4] = { {-1000.0,-1000.0,-1000.0,-1000.0},{ -1000.0,-1000.0,-1000.0,-1000.0 } ,{ -1000.0,-1000.0,-1000.0,-1000.0 } ,{ -1000.0,-1000.0,-1000.0,-1000.0 } ,{ -1000.0,-1000.0,-1000.0,-1000.0 } };
		const Configuration configs[nThreads] = { conf.at(0),conf.at(1),conf.at(2),conf.at(3),conf.at(4) };
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
#pragma omp parallel num_threads(nThreads) shared(buffData, xcoords,configs,loggersvector)
				{

					int threadid = omp_get_thread_num();
					//std::vector <double> xcoord = xcoords[threadid];
					Configuration config = configs[threadid];

					const LoggerParameters lP = config.loggerParameters;
					const ModelParameters mP = config.modelParameters;
					InitialConditions iC = config.initialConditions;

					// configurate force object
					PotentialForce potentialForce;
					potentialForce.E = E;
					potentialForce.G = mP.G;
					potentialForce.L = mP.L;
					potentialForce.powsigma = pow(mP.sigma, 2.0);
					double xMT;
					double xBeadl;
					double xBeadr;
					double xMol;

					if ((xcoords[threadid][0] == -1000.0) && (xcoords[threadid][1] == -1000.0) && (xcoords[threadid][2] == -1000.0) && (xcoords[threadid][3] == -1000.0)) {
						xMT = iC.xMT;
						xBeadl = iC.xBeadl;
						xBeadr = iC.xBeadr;
						xMol = iC.xMol;
					}
					else {
						xMT = xcoords[threadid][0];
						xBeadl = xcoords[threadid][1];
						xBeadr = xcoords[threadid][2];
						xMol = xcoords[threadid][3];
					}
					for (int iter = 0; iter < 900'000/3; iter ++) {

						const double MT_Mol_force = potentialForce.calc(xMol - xMT);
						//xMT = -0.5*((force( xMT-xMol) / MTstiff) - xBeadl - xBeadr);
						//const double next_xMT =   ((MTstiffL*xBeadl) + (MTstiffR*xBeadr)-force(xMT - xMol))/(MTstiffL+ MTstiffR);
						const double next_xMT = xMT + (sP.expTime / mP.gammaMT)*(((-mP.MTstiffL)*(xMT - xBeadl)) + (mP.MTstiffR*(xBeadr - xMT)) - (MT_Mol_force));
						const double next_xBeadl = xBeadl + (sP.expTime / mP.gammaBead)*(((-mP.trapstiff)*(xBeadl - iC.xTrapl)) + (mP.MTstiffL*(xMT - xBeadl))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter]);

						const double next_xBeadr = xBeadr + (sP.expTime / mP.gammaBead)*(((-mP.MTstiffR)*(xBeadr - xMT)) + ((-mP.trapstiff)*(xBeadr - iC.xTrapr))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter + (900'000/3)]);

						const double next_xMol = xMol + (sP.expTime / mP.gammaMol) *(MT_Mol_force + mP.molstiff*(iC.xPed - xMol)) + sqrt(2.0*mP.DMol*sP.expTime) *(buffData[iter + (2*900'000 / 3)]);

						xMT = next_xMT;
						xBeadl = next_xBeadl;
						xBeadr = next_xBeadr;
						xMol = next_xMol;
					}
					/*xcoord.push_back(xMT);
					xcoord.push_back(xBeadl);
					xcoord.push_back(xBeadr);
					xcoord.push_back(xMol);
					xcoords[threadid] = xcoord;
					*/
					//int threadid = omp_get_thread_num();
					xcoords[threadid][0] = xMT;
					xcoords[threadid][1] = xBeadl;
					xcoords[threadid][2] = xBeadr;
					xcoords[threadid][3] = xMol;
					loggersvector.at(threadid).save(xcoords[threadid][0], xcoords[threadid][1], xcoords[threadid][2], xcoords[threadid][3]);
				}
			}
			/*
			for (int it = 0; it <= 4; it++) {
				loggersvector.at(it).save(xcoords[it][0], xcoords[it][1], xcoords[it][2], xcoords[it][3]);
			}
			*/
			//std::cout << xcoords[0][1] << std::endl;
			if (savedstep % 100 == 0) {
				double procent = round(100*100 * savedstep / (100'000))/100;
				std::cout << procent << "%" << std::endl;
				std::cout << __rdtsc() << std::endl;
				//std::cout << nst << std::endl;
			}
		}
		for (int it = 0; it <= 4; it++) {
			loggersvector.at(it).closeLogger();
		}
		
		const auto deletionResult = vslDeleteStream(&stream);
		if (deletionResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "error deleting stream" };
		}

		
		//////////////////////
		////////////////////
		//////////////////////
	}
	if (mode == 2) {
			if (cmdOptionExists(argv, argv + argc, "-h"))
			{
				// Do stuff
				std::cout << "Sorry users, no help donations today." << std::endl;
			}
			char * param_input_filename = getCmdOption(argv, argv + argc, "-paramsfile");
			char * output_filename = getCmdOption(argv, argv + argc, "-resultfile");

			//////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////
			////////////////// Random Number Generator  ////////////////

			VSLStreamStatePtr stream = nullptr;
			const auto creationResult = vslNewStream(&stream, VSL_BRNG_SFMT19937, 777);
			if (creationResult != VSL_STATUS_OK) {
				throw std::runtime_error{ "can't create VSL stream" };
			}

			static constexpr std::size_t buffSize = 900000;// Maximum size should fit L3 chache //1200'000
			int intbufferSize = 900000;
			//std::size_t buffSize = 1200'000;
			std::vector<double> buff(buffSize);

			const auto generateNumbers = [&buff, stream]() {

				static constexpr std::size_t nPerThread = buffSize / nThreads;
				static_assert(buffSize == nPerThread * nThreads, "buffsize must be multiple of nThreads");

				const auto buffData = buff.data();
				const auto _stream = stream;
#pragma omp parallel num_threads(nThreads) shared(buffData, _stream)
				{
					const std::size_t begin = nPerThread * omp_get_thread_num();
					const auto generationResult = vdRngGaussian(
						VSL_RNG_METHOD_GAUSSIAN_ICDF, _stream,
						nPerThread, buffData + begin,
						0.0, 1.0
					);
					if (generationResult != VSL_STATUS_OK) {
						throw std::runtime_error{ "can't generate numbers" };
					}
				}

			};

			////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////












			//const ModelParameters mP;
			//InitialConditions iC;
			//const SimulationParameters sP;

			SimulationParameters sP = load_simulationparams(inputparamfile);
			std::vector <Configuration> conf = load_configuration(inputparamfile);

			//const auto logger = createLogger(LoggerType::BINARY_FILE, lP);

			BinaryFileLogger logger0 = create_logger(conf.at(0).loggerParameters);
			BinaryFileLogger logger1 = create_logger(conf.at(1).loggerParameters);
			BinaryFileLogger logger2 = create_logger(conf.at(2).loggerParameters);
			BinaryFileLogger logger3 = create_logger(conf.at(3).loggerParameters);
			BinaryFileLogger logger4 = create_logger(conf.at(4).loggerParameters);

			///////////////////////////////
			//////////////////////////////
			///////////// Iterations
			///////////////////////////
			/*
			const int frequency = (int)round(sP.microsteps);
			const int totalsimsteps = (int)round((sP.microsteps)*(sP.nTotal));
			const int nsteps = (int)round(3*(totalsimsteps / frequency) / intbufferSize);
			*/

			double xcoords[nThreads][4] = { { -1000.0,-1000.0,-1000.0,-1000.0 },{ -1000.0,-1000.0,-1000.0,-1000.0 },{ -1000.0,-1000.0,-1000.0,-1000.0 },{ -1000.0,-1000.0,-1000.0,-1000.0 },{ -1000.0,-1000.0,-1000.0,-1000.0 } };
			const Configuration configs[nThreads] = { conf.at(0),conf.at(1),conf.at(2),conf.at(3),conf.at(4) };
			/*
			int frequency = sP.microsteps;
			int totalsimsteps = sP.microsteps*(sP.nTotal);
			int nsteps = 3 * (totalsimsteps / frequency) / intbufferSize;
			int dd = 2;

			frequency 100'000
			totalsimsteps 100'000'000'000
			buffer 1'000'000
			randms per simstep 3
			*/
			int count = 0;
			for (int macrostep = 0; macrostep <  (900'000 / 900'000); macrostep++) {
					generateNumbers();
					const auto buffData = buff.data();
					
					for(int threadid=0; threadid<nThreads; threadid++){
						//int threadid = omp_get_thread_num();
						//std::vector <double> xcoord = xcoords[threadid];
						Configuration config = configs[threadid];

						const LoggerParameters lP = config.loggerParameters;
						const ModelParameters mP = config.modelParameters;
						InitialConditions iC = config.initialConditions;

						// configurate force object
						PotentialForce potentialForce;
						potentialForce.E = E;
						potentialForce.G = mP.G;
						potentialForce.L = mP.L;
						potentialForce.powsigma = pow(mP.sigma, 2.0);
						double xMT;
						double xBeadl;
						double xBeadr;
						double xMol;
						
						if ((xcoords[threadid][0] == -1000.0) && (xcoords[threadid][1] == -1000.0) && (xcoords[threadid][2] == -1000.0) && (xcoords[threadid][3] == -1000.0)) {
							xMT = iC.xMT;
							xBeadl = iC.xBeadl;
							xBeadr = iC.xBeadr;
							xMol = iC.xMol;
							count++;
						}
						else {
							xMT = xcoords[threadid][0];
							xBeadl = xcoords[threadid][1];
							xBeadr = xcoords[threadid][2];
							xMol = xcoords[threadid][3];
						}
						for (int iter = 0; iter < 300'000; iter++ ) {

							const double MT_Mol_force = potentialForce.calc(xMol - xMT);
							//xMT = -0.5*((force( xMT-xMol) / MTstiff) - xBeadl - xBeadr);
							//const double next_xMT =   ((MTstiffL*xBeadl) + (MTstiffR*xBeadr)-force(xMT - xMol))/(MTstiffL+ MTstiffR);
							const double next_xMT = xMT + (sP.expTime / mP.gammaMT)*(((-mP.MTstiffL)*(xMT - xBeadl)) + (mP.MTstiffR*(xBeadr - xMT)) - (MT_Mol_force));
							const double next_xBeadl = xBeadl + (sP.expTime / mP.gammaBead)*(((-mP.trapstiff)*(xBeadl - iC.xTrapl)) + (mP.MTstiffL*(xMT - xBeadl))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter]);

							const double next_xBeadr = xBeadr + (sP.expTime / mP.gammaBead)*(((-mP.MTstiffR)*(xBeadr - xMT)) + ((-mP.trapstiff)*(xBeadr - iC.xTrapr))) + sqrt(2.0*mP.DBead*sP.expTime)*(buffData[iter + 300'000]);

							const double next_xMol = xMol + (sP.expTime / mP.gammaMol) *(MT_Mol_force + mP.molstiff*(iC.xPed - xMol)) + sqrt(2.0*mP.DMol*sP.expTime) *(buffData[iter +600'000]);

							xMT = next_xMT;
							xBeadl = next_xBeadl;
							xBeadr = next_xBeadr;
							xMol = next_xMol;
						
							xcoords[threadid][0] = xMT;
							xcoords[threadid][1] = xBeadl;
							xcoords[threadid][2] = xBeadr;
							xcoords[threadid][3] = xMol;

							if(threadid==0) {
								logger0.save(xcoords[0][0], xcoords[0][1], xcoords[0][2], xcoords[0][3]);
								
						}
							if (threadid == 1) {
								logger1.save(xcoords[1][0], xcoords[1][1], xcoords[1][2], xcoords[1][3]);

							}
							if (threadid == 2) {
								logger2.save(xcoords[2][0], xcoords[2][1], xcoords[2][2], xcoords[2][3]);

							}
							if (threadid == 3) {
								logger3.save(xcoords[3][0], xcoords[3][1], xcoords[3][2], xcoords[3][3]);

							}
							if (threadid == 4) {
								logger4.save(xcoords[4][0], xcoords[4][1], xcoords[4][2], xcoords[4][3]);

							}
							
							
							
						
						
						}
						

						
						
					}
					std::cout << "macrostep" << macrostep << "count" << count << std::endl;
			}

				

				
				
			
			const auto deletionResult = vslDeleteStream(&stream);
			if (deletionResult != VSL_STATUS_OK) {
				throw std::runtime_error{ "error deleting stream" };
			}


			//////////////////////
			////////////////////
			//////////////////////
			logger0.closeLogger();
			logger1.closeLogger();
			logger2.closeLogger();
			logger3.closeLogger();
			logger4.closeLogger();
		}
return 0;
}

//std::chrono
//_rdtscp  in #immintrin.h