#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <array>
#include <cstdint>
#include "emissionprobabilitycomputer.hpp"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(shared_ptr<UniqueKmers> uniquekmers, ProbabilityTable* probabilities)
	:uniquekmers(uniquekmers),
	 probabilities(probabilities),
	 all_zeros(true)
{
	vector<unsigned short> unique_alleles;
	uniquekmers->get_allele_ids(unique_alleles);
	if (unique_alleles.empty()) return;
	unsigned short max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
	this->state_to_prob = vector< vector<double>>(max_allele+1, vector<double>(max_allele+1));

	// stage5 optimization: precompute, ONCE per column, two tables that the inner
	// (a1, a2) x M loop would otherwise re-derive on every iteration:
	//   P_table[i][ec]    = probability of read-kmer-count for kmer i under expected_kmer_count ec in {0,1,2}
	//   kmer_bit[i][a]    = whether kmer i is present on allele a
	// This collapses the hot multiplication chain from
	//   probabilities->get_probability(cov, readcount).get_probability_of(ec)   (~30 cycles, 2x .at())
	// down to a single contiguous-array lookup (~3 cycles).
	const size_t M = uniquekmers->size();
	const unsigned short cov = uniquekmers->get_coverage();

	vector<array<double, 3>> P_table(M);
	for (size_t i = 0; i < M; ++i) {
		CopyNumber cn = probabilities->get_probability(cov, uniquekmers->get_readcount_of(i));
		P_table[i][0] = cn.get_probability_of(0);
		P_table[i][1] = cn.get_probability_of(1);
		P_table[i][2] = cn.get_probability_of(2);
	}

	// kmer_bit is M x (max_allele+1) packed into a flat vector for cache-friendly access.
	const size_t row_stride = (size_t)max_allele + 1;
	vector<uint8_t> kmer_bit(M * row_stride, 0);
	for (size_t i = 0; i < M; ++i) {
		for (auto a : unique_alleles) {
			kmer_bit[i * row_stride + a] = uniquekmers->kmer_on_allele(i, a) ? 1 : 0;
		}
	}

	// Hoist undefined-allele lookups outside the M loop.
	vector<bool> undef(max_allele + 1, false);
	for (auto a : unique_alleles) {
		undef[a] = uniquekmers->is_undefined_allele(a);
	}

	const double one_third = 1.0 / 3.0;
	const double half = 0.5;

	for (auto a1 : unique_alleles) {
		const bool a1_und = undef[a1];
		for (auto a2 : unique_alleles) {
			const bool a2_und = undef[a2];
			double result = 1.0;
			if (a1_und && a2_und) {
				// All three copy-number bins equally weighted.
				for (size_t i = 0; i < M; ++i) {
					result *= one_third * (P_table[i][0] + P_table[i][1] + P_table[i][2]);
				}
			} else if (a1_und || a2_und) {
				// One side undefined: average between expected_kmer_count and expected_kmer_count+1.
				// expected_kmer_count must be < 2 (asserted upstream), so ec+1 <= 2 in P_table[].
				for (size_t i = 0; i < M; ++i) {
					const int ec = (int)kmer_bit[i * row_stride + a1] + (int)kmer_bit[i * row_stride + a2];
					assert(ec < 2);
					result *= half * (P_table[i][ec] + P_table[i][ec + 1]);
				}
			} else {
				// Defined alleles: expected_kmer_count is exact.
				for (size_t i = 0; i < M; ++i) {
					const int ec = (int)kmer_bit[i * row_stride + a1] + (int)kmer_bit[i * row_stride + a2];
					result *= P_table[i][ec];
				}
			}
			this->state_to_prob[a1][a2] = result;
			if (result > 0) this->all_zeros = false;
		}
	}

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

double EmissionProbabilityComputer::get_emission_probability(unsigned short allele_id1, unsigned short allele_id2) const {
	if (this->all_zeros) return 1.0;
	return this->state_to_prob[allele_id1][allele_id2];
}

double EmissionProbabilityComputer::compute_emission_probability(unsigned short allele_id1, unsigned short allele_id2, bool a1_undefined, bool a2_undefined){
	double result = 1.0;
	for (size_t i = 0; i < this->uniquekmers->size(); ++i){
		unsigned int expected_kmer_count = this->uniquekmers->kmer_on_allele(i, allele_id1) + this->uniquekmers->kmer_on_allele(i, allele_id2);
		if (a1_undefined && a2_undefined) {
			// all kmers can have copy numbers 0-2
			result *= (1.0 / 3.0) * (this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(0) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(1) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(2));
		} else if (a1_undefined || a2_undefined) {
			// two possible copy numbers
			assert (expected_kmer_count < 2);
			result *= 0.5 * (this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count + 1));
		} else {
			// expected kmer count is known
			result *= this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count);
		}
	}
	return result;
}
