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
	double kOn2, kOff2;

	double** transitionMatrix1{ assign_rates(numStates, kOn1, kOff1) };
	double** transitionMatrix2{ assign_rates(numStates, kOn2, kOff2) };

	bool bindingDynamics = true;
	double B = 0.0;

	static double** assign_rates(int numStates, double kOn, double kOff) {

		double** transitionMatrix = new double*[numStates];
		for (int i = 0; i < numStates; ++i) {
			transitionMatrix[i] = new double[numStates];
		}

		transitionMatrix[0][0] = -kOn;
		transitionMatrix[0][1] = kOn;

		transitionMatrix[1][1] = kOff;
		transitionMatrix[1][2] = -kOff;

		return transitionMatrix;
	}
};

class MoleculeState
{
public:
	unsigned short label = 0;
	double xMol = 0.0;
	double logpotentialForce = 0.0;
	double binding = 0.0; // binding of the second MOLECULE, first mol is always bound
	//double currentWell = 0.0;
	double phi = 0.0;
	double potTorque = 0.0;
	double deltaG = 0.0;
	double MToffset = 0.0;

	std::vector<double> expRands;
	std::vector<double> livingTimes;

	MoleculeState(unsigned short lbl, double mt_offset) : label(lbl), MToffset(mt_offset) {	};

	template <typename F>
	void iterateFields(F&& f) {
		std::string strLabel = std::to_string(label);

		f(&MoleculeState::xMol, "xMol" + strLabel);
		f(&MoleculeState::logpotentialForce, "logpotentialForce" + strLabel);
		f(&MoleculeState::binding, "binding" + strLabel);
		//f(&MoleculeState::currentWell, "currentWell" + strLabel);
		f(&MoleculeState::phi, "phi" + strLabel);
		f(&MoleculeState::potTorque, "potTorque" + strLabel);
		f(&MoleculeState::deltaG, "deltaG" + strLabel);
	}
};

struct SystemState
{
	double xMT;
	double xBeadl;
	double xBeadr;
	double xTrapl; 
	double xTrapr; 
	double Time = 0.0;
	double direction = 1.0;

	MoleculeState firstMol{ 1, 0 };
	MoleculeState secondMol{ 2, 0.0006 };

	//#pragma omp declare simd
	template <typename F>
	static void iterateFields(F&& f) {
		f(&SystemState::xMT, "xMT");
		f(&SystemState::xBeadl, "xBeadl");
		f(&SystemState::xBeadr, "xBeadr");

		f(&SystemState::xTrapl, "xTrapl");
		f(&SystemState::xTrapr, "xTrapr");
		f(&SystemState::Time, "Time");
		f(&SystemState::direction, "direction");
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
