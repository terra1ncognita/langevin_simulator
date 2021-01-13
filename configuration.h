#pragma once
#include <string>
#include <omp.h>
#include <vector>
#include <string>

// Global parameters of simulations -> define iterations and steps
struct SimulationParameters
{
	double expTime = 1e-12;
	double simulationTime = 3e-3;
	unsigned int iterationsbetweenSavings = 15'000'000;
	unsigned int iterationsbetweenTrapsUpdate = 15'000'000;
	unsigned int totalsavings = int((simulationTime / expTime) / iterationsbetweenSavings);
	
	unsigned int randomsPeriter = 4;
	std::size_t buffsize = 15'000'000 * randomsPeriter;
	unsigned int stepsperbuffer = static_cast<unsigned int>(std::floor(buffsize / randomsPeriter));
	
	unsigned int savingsPerMacrostep = stepsperbuffer / iterationsbetweenSavings;
	unsigned int macrostepMax = totalsavings / savingsPerMacrostep;
	unsigned int trapsUpdateTest = iterationsbetweenTrapsUpdate / iterationsbetweenSavings;

	unsigned short rndThreads = 4;
};


// Classes for objects that store simulation configs: LoggerParameters, ModelParameters, InitialConditions
struct LoggerParameters
{
	//enum class FilenameTemplate { PrefixName };
	//FilenameTemplate filenametemplate;
	std::string filepath;
	std::string name;
};

struct ModelParameters
{
	/////
	//Global paramters
	std::string name = "";
	double T;					//temperature
	double kT;

	//Parameters of potential
	double G, G2;					// (* kT | Depth of the potential *)
	double L;					//(* um | period of the periodic potential *)
	double sigma;			//(* um | width of the binding well *)
	double A;           // width of asymmetric potential, um
	double m;           // center of well, um
									//Parameters of diffusion
	double DMol;					//(* um^2/s | free diffusion coefficient of the protein in water *)
	double DBeadL;					// (* um^2/s | free diffusion coefficient of 0.5 um bead in water *)
	double DBeadR;
	double DMT;

	double gammaMol;		//(* pN s/um | friction drag coefficient for protein *)
	double gammaBeadL;		//(* pN s/um | friction drag coefficient for 0.5 um bead *)
	double gammaBeadR;
	double gammaMT;				//(* pN s/um | friction drag coefficient for 0.5 um for MT *)

	double gammaQuasiviscous;
											// Parameters of stiffness
	double trapstiffL;			//(* pN/um | stiffness of the trap *) 
	double trapstiffR;

	double MTstiffWeakSlopeL;
	double MTstiffWeakBoundaryL;
	double MTstiffParabolicAL;
	double MTstiffParabolicBL;
	double MTstiffParabolicCL;
	double MTstiffStrongBoundaryL;
	double MTstiffStrongSlopeL;
	double MTstiffStrongIntersectL;

	double MTstiffWeakSlopeR;
	double MTstiffWeakBoundaryR;
	double MTstiffParabolicAR;
	double MTstiffParabolicBR;
	double MTstiffParabolicCR;
	double MTstiffStrongBoundaryR;
	double MTstiffStrongSlopeR;
	double MTstiffStrongIntersectR;

	double molStiffWeakSlope;
	double molStiffBoundary;
	double molStiffStrongSlope;

	double MTlength;
	double molstiff;				//(*pN / um| stiffness of the NDC80 *)
	double feedbackFreq;
	double DmblMoveAmplitude;
	double prestretchTotalForce;
	double movementTotalForce;

	double rotFriction;
	double rotStiffness;
	double molLength;
	double domainsDistance;
	double rotWellWidth;
	double rotWellDepth;
	double iniPhi;

	int numStates = 2;

	double kOn1, kOff1;
	// double kOn2, kOff2;

	double** transitionMatrix;

	bool bindingDynamics = true;
	double B = 0.0;
};

struct SystemState
{
	double xMol;
	double xMT;
	double xBeadl;
	double xBeadr;
	double xTrapl; 
	double xTrapr; 
	double Time = 0.0;
	double direction = 1.0;
	double logpotentialForce;

	double binding = 0.0;
	double currentWell = 0.0;

	double phi = 0.0;
	double potTorque = 0.0;
	double deltaG = 0.0;

	//#pragma omp declare simd
	template <typename F>
	static void iterateFields(F&& f) {
		f(&SystemState::xMol, "xMol");
		f(&SystemState::xMT, "xMT");
		f(&SystemState::xBeadl, "xBeadl");
		f(&SystemState::xBeadr, "xBeadr");

		f(&SystemState::xTrapl, "xTrapl");
		f(&SystemState::xTrapr, "xTrapr");
		f(&SystemState::Time, "Time");
		f(&SystemState::direction, "direction");
		f(&SystemState::logpotentialForce, "logpotentialForce");

		f(&SystemState::binding, "binding");
		f(&SystemState::phi, "phi");
		f(&SystemState::potTorque, "potTorque");
		f(&SystemState::deltaG, "deltaG");
	}
};

struct InitialConditions
{
	SystemState initialState;
};

// Composition of parameters
struct Configuration
{
	SimulationParameters simulationParameters;
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
