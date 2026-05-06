#include <cmath>
#include <cassert>
#include "samplingtransitions.hpp"

SamplingTransitions::SamplingTransitions(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, double effective_N) {
	assert(from_variant <= to_variant);
	// using same formula as in WhatsHap
//	double distance = (to_variant - from_variant) * 0.000001 * ((double) recomb_rate) * 4.0 * effective_N;
	double distance = (to_variant - from_variant) * 0.000004 * ((double) recomb_rate) * effective_N;
	// use Li-Stephans pair HMM transitions TODO: correct?
	double recomb_prob = (1.0 - exp(-distance / (double) nr_paths) )* (1.0 / (double) nr_paths);
	this->cost = -10.0 * log10(recomb_prob);

//	this->cost = 10.0;
}

unsigned int SamplingTransitions::compute_transition_cost(bool recombination) {
	if (recombination) {
		return this->cost;
	} else {
		return 0;
	}
}
