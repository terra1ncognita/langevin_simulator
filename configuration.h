#pragma once
#include <string>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>
#include "json.hpp"
#include "polynomial.h"

using json = nlohmann::json;

class ForceExtensionCurve {
public:

	ForceExtensionCurve() {};

	ForceExtensionCurve(LaurentPolynomial left, LaurentPolynomial right, double inflection_point) 
		: left{ left }, right{ right }, inflection_point{inflection_point} 
	{ check_validity(); }

	ForceExtensionCurve(const json& jsObj, std::string key) {
		check_key_json(jsObj, key);

		left = LaurentPolynomial(jsObj[key], "left");
		right = LaurentPolynomial(jsObj[key], "right");
		inflection_point = static_cast<double>(jsObj[key]["inflection_point"]);

		check_validity();
	}

	double operator() (double y) const {
		double x = y - inflection_point;
		if (x < left.pole() || x > right.pole()) {
			std::stringstream ss;
			ss << "Displacement " << x << " is out of boundaries [" << left.pole() << ", " << right.pole() << "]";
			throw std::runtime_error(ss.str());
		}
		if (x < 0) {
			return left(x);
		}
		return right(x);
	}

	double left_pole() const {
		return left.pole();
	}

	double right_pole() const {
		return right.pole();
	}

private:
	double inflection_point;
	LaurentPolynomial left, right;

	void check_validity() {
		if (left.pole() >= right.pole()) {
			throw std::runtime_error("Left pole is greater than right pole");
		}
		if (inflection_point >= right.pole() || inflection_point <= left.pole()) {
			throw std::runtime_error("Inflection point is not between poles");
		}
	}
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

	ForceExtensionCurve fec;

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
