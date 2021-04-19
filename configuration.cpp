#include "configuration.h"
#include <iostream>


std::ostream& operator<< (std::ostream &out, const SimulationParameters& sp) {
	out << "SimulationParameters(";
	out << "expTime=" << sp.expTime << ", ";
	out << "simulationTime=" << sp.simulationTime << ", ";
	out << "iterationsbetweenSavings=" << sp.iterationsbetweenSavings << ", ";
	out << "iterationsbetweenTrapsUpdate=" << sp.iterationsbetweenTrapsUpdate << ", ";
	out << "totalsavings=" << sp.totalsavings << ", ";

	out << "randomsPeriter=" << sp.randomsPeriter << ", ";
	out << "buffsize=" << sp.buffsize << ", ";
	out << "stepsperbuffer=" << sp.stepsperbuffer << ", ";

	out << "savingsPerMacrostep=" << sp.savingsPerMacrostep << ", ";
	out << "macrostepMax=" << sp.macrostepMax << ", ";
	out << "trapsUpdateTest=" << sp.trapsUpdateTest << ", ";

	out << "freeMotionTime=" << sp.freeMotionTime << ", ";
	out << "macrostepsFree=" << sp.macrostepsFree << ", ";

	out << "rndThreads=" << sp.rndThreads;
	out << ")";
	return out;
}
