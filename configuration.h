#pragma once
#include <string>
# include <omp.h>
#include <vector>
#include <string>

// Global parameters of simulations -> define iterations and steps
struct SimulationParameters
{
	// simulation parameters
	double expTime;//
	int microsteps;// 10MHz scanning
	int nTotal ;	//10^-2 seconds total
};


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
	//Global paramters
	double T ;					//temperature
	double kT;

	//Parameters of potential
	double G ;					// (* kT | Depth of the potential *)
	double L ;					//(* um | period of the periodic potential *)
	double sigma;			//(* um | width of the binding well *)

									//Parameters of diffusion
	double DMol ;					//(* um^2/s | free diffusion coefficient of the protein in water *)
	double DBead ;					// (* um^2/s | free diffusion coefficient of 0.5 um bead in water *)
	double DMT;

	double gammaMol;		//(* pN s/um | friction drag coefficient for protein *)
	double gammaBead ;		//(* pN s/um | friction drag coefficient for 0.5 um bead *)
	double gammaMT ;				//(* pN s/um | friction drag coefficient for 0.5 um for MT *)

											// Parameters of stiffness
	double trapstiff;			//(* pN/um | stiffness of the trap *) 
	double MTstiffL ;			//(* pN/um | stiffness of the MT *) 
	double MTstiffR;
	double molstiff ;				//(*pN / um| stiffness of the NDC80 *)

};

struct SystemState
{
	double xMol;
	double xMT;
	double xBeadl;
	double xBeadr;

	template <typename F>
	static void iterateFields(F&& f) {
		f(&SystemState::xMol, "xMol");
		f(&SystemState::xMT, "xMT");
		f(&SystemState::xBeadl, "xBeadl");
		f(&SystemState::xBeadr, "xBeadr");
	}
};

struct InitialConditions
{
	SystemState initialState;

	double xPed;   ////////////// Is it really iC????
	double xTrapl; // Must be negative for prestretch ////////////// Is it really iC????
	double xTrapr; // Must be positive for prestretch ////////////// Is it really iC????
};
// Composition of parameters
struct Configuration
{
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
