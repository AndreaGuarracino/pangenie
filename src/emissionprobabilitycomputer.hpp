#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "uniquekmers.hpp"
#include "copynumber.hpp"
#include "columnindexer.hpp"
#include "probabilitytable.hpp"

/** 
* Computes the emission probabilities for a variant position.
**/

typedef std::vector<std::vector<double>> ProbabilityMatrix;

class EmissionProbabilityComputer {
public:
	/**
	* @param uniquekmers all unique kmers for this position
	 **/
	EmissionProbabilityComputer(const std::shared_ptr<UniqueKmers>& uniquekmers, ProbabilityTable* probabilities);
	/** get emission probability for a state in the HMM **/
	double get_emission_probability(unsigned short allele_id1, unsigned short allele_id2) const noexcept {
		if (this->all_zeros) return 1.0;
		return this->state_to_prob[(size_t)allele_id1 * ((size_t)this->max_allele + 1) + allele_id2];
	}

private:
	UniqueKmers* uniquekmers;
	ProbabilityTable* probabilities;
	bool all_zeros;
	unsigned short max_allele;
        std::vector<double> state_to_prob;
	double compute_emission_probability(unsigned short allele1, unsigned short allele2, bool allele1_undefined, bool allele2_undefined);
};
# endif // EMISSIONPROBABILITYCOMPUTER_H
