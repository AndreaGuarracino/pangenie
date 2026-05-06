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
	double recomb_prob = (1.0 - exp(-distance / (double) nr_paths) )* (1.0 / (double) nr_paths);
	double no_recomb_prob = exp(-distance / (double) nr_paths) + recomb_prob;
	this->probabilities = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};
}

double TransitionProbabilityComputer::compute_transition_prob(unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4){
	if (this->uniform) {
		return 1.0;
	} else {
		// determine number of recombination events
		unsigned int nr_switches = 0;
		if (path_id1 != path_id3) nr_switches += 1;
		if (path_id2 != path_id4) nr_switches += 1;
		return this->probabilities[nr_switches];
	}
}

double TransitionProbabilityComputer::compute_transition_prob(unsigned short nr_switches) {
	if (this->uniform) {
		return 1.0;
	} else {
		return this->probabilities[nr_switches];
	}
}
