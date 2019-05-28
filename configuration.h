#pragma once
#include <string>
# include <omp.h>
#include <vector>
#include <string>

// Global parameters of simulations -> define iterations and steps
//struct SimulationParameters
//{
	// simulation parameters
//	double expTime;//
	//int microsteps;// 10MHz scanning
	//int nTotal ;	//10^-2 seconds total
//};


// Classes for objects that store simulation configs: LoggerParameters, ModelParameters, InitialConditions
struct LoggerParameters
{
	enum class FilenameTemplate { PrefixName };
	FilenameTemplate filenametemplate;
	std::string filepath;
	std::string name;
};
struct ModelParameters
{
	/////
	double expTime;
	//Global paramters
	double T ;					//temperature
	double kT;

	//Parameters of potential
	double G ;					// (* kT | Depth of the potential *)
	double L ;					//(* um | period of the periodic potential *)
	double sigma;			//(* um | width of the binding well *)
	double A;           // width of asymmetric potential, um
	double m;           // center of well, um
									//Parameters of diffusion
	double DMol ;					//(* um^2/s | free diffusion coefficient of the protein in water *)
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
	double MTstiffL ;			//(* pN/um | stiffness of the MT *) 
	double MTstiffR;
	double MTlowstiff;
	double MTrelaxedLengthL;
	double MTrelaxedLengthR;

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
	double molstiff ;				//(*pN / um| stiffness of the NDC80 *)
	double feedbackFreq;
	double DmblMoveAmplitude;
	double prestretchTotalForce;
	double movementTotalForce;
};

struct SystemState
{
	double xMol;
	double xMT;
	double xBeadl;
	double xBeadr;
	double xTrapl; 
	double xTrapr; 
	double Time=0.0;
	double direction = 1.0;
	double logpotentialForce;

	double vMol = 0.0;
	double vMT = 0.0;

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
	}
};

struct InitialConditions
{
	SystemState initialState;
	
	double xPed;   ////////////// Is it really iC????
	//double xTrapl; // Must be negative for prestretch ////////////// Is it really iC????
	//double xTrapr; // Must be positive for prestretch ////////////// Is it really iC????
};
// Composition of parameters
struct Configuration
{
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
