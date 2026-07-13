#include <cassert>
#include <math.h>
#include "transitionprobabilitycomputer.hpp"
#include <iostream>

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, bool uniform, double effective_N)
{
	assert(from_variant <= to_variant );
	this->uniform = uniform;
	// using same formula as in WhatsHap
//	double distance = (to_variant - from_variant) * 0.000001 * ((double) recomb_rate) * 4.0 * effective_N;
	double distance = (to_variant - from_variant) * 0.000004 * ((double) recomb_rate) * effective_N;
	// use Li-Stephans pair HMM transitions TODO: correct?
	double no_switch = exp(-distance / (double) nr_paths);
	double recomb_prob = (1.0 - no_switch) * (1.0 / (double) nr_paths);
	double no_recomb_prob = no_switch + recomb_prob;
	this->probabilities = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};
}
