#ifndef TRANSITIONPROBABILITYCOMPUTER_HPP
#define TRANSITIONPROBABILITYCOMPUTER_HPP

#include <array>
#include "uniquekmers.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class TransitionProbabilityComputer {
public:
	TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, bool uniform = false, double effective_N = 25000.0);
	// computes the transition probability given previous and current paths
	double compute_transition_prob(unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4) const noexcept {
		if (uniform) return 1.0;
		const unsigned int nr_switches = (path_id1 != path_id3) + (path_id2 != path_id4);
		return probabilities[nr_switches];
	}
	// computes the transition probability based on the number of haplotype switches
	double compute_transition_prob(unsigned short nr_switches) const noexcept {
		return uniform ? 1.0 : probabilities[nr_switches];
	}
private:
	std::array<double, 3> probabilities;
	bool uniform;
	
};
#endif // TRANSITIONPROBABILITYCOMPUTER_HPP
