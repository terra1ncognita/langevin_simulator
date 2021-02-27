#pragma once
#include <string>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>

class LaurentPolynomial {
public:
	LaurentPolynomial() : n{ 0 }, m{ 0 }, pole{ 0 }, coeffs{}  {}
	LaurentPolynomial(int n, int m, double pole, std::vector<double> coeffs) : n(n), m(m), pole(pole), coeffs(coeffs) {
		if (n - m + 1 != coeffs.size()) {
			throw std::runtime_error("n, m and coeff size are not in agreement!");
		}
	}

	double operator()(double x) const {
		double res = coeffs[-m];
		for (int i = 0; i < - m; ++i) {
			res += coeffs[i] * pow(x - pole, i + m);
		}
		for (int i = 1-m; i < n-m+1; ++i) {
			res += coeffs[i] * pow(x, i + m);
		}
		return res;
	}
private:
	int n, m;
	double pole;
	std::vector<double> coeffs;
};

class RationalFunction {
public:
	RationalFunction() : p{}, q{} {}
	RationalFunction(std::vector<double> pp, std::vector<double> qq) : p{ pp }, q{ qq } {}

	static double evaluate_polynomial(const std::vector<double> &coeffs, double x) {
		double curr_pow = 1.0, res = 0.0;
		for (int i = 0; i < coeffs.size(); ++i) {
			res += curr_pow * coeffs[i];
			curr_pow *= x;
		}
		return res;
	}

	double operator() (double x) const {
		double num = evaluate_polynomial(p, x);
		double denom = evaluate_polynomial(q, x);
		return num / denom;
	}
private:
	std::vector<double> p, q;
};


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

	double freeMotionTime = 5e-3;
	unsigned int macrostepsFree = static_cast<unsigned int>(ceil(freeMotionTime / expTime / iterationsbetweenSavings / savingsPerMacrostep));

	unsigned short rndThreads = 4;
};


std::ostream& operator<< (std::ostream &out, const SimulationParameters& sp);


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

	LaurentPolynomial MTextension;
	RationalFunction MTcompression;

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
	double xMT;
	double xBeadl;
	double xBeadr;
	double xTrapl; 
	double xTrapr; 
	double Time = 0.0;
	double direction = 1.0;

	bool instaBindingDynamics = false;
	bool isFree = true;
	double xMTbinding;

	double binding = 0.0;
	double currentWell = 0.0;

	
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

		f(&SystemState::binding, "binding");
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
