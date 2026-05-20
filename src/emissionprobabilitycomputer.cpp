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
        this->max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
        size_t stride = (size_t)this->max_allele + 1;
        this->state_to_prob.assign(stride * stride, 0.0);

        const size_t M = uniquekmers->size();
        const unsigned short cov = uniquekmers->get_coverage();

        vector<array<double, 3>> P_table(M);
        for (size_t i = 0; i < M; ++i) {
                CopyNumber cn = probabilities->get_probability(cov, uniquekmers->get_readcount_of(i));
                P_table[i][0] = cn.get_probability_of(0);
                P_table[i][1] = cn.get_probability_of(1);
                P_table[i][2] = cn.get_probability_of(2);
        }

        const size_t row_stride = stride;
        vector<uint8_t> kmer_bit(M * row_stride, 0);
        for (size_t i = 0; i < M; ++i) {
                for (auto a : unique_alleles) {
                        kmer_bit[i * row_stride + a] = uniquekmers->kmer_on_allele(i, a) ? 1 : 0;
                }
        }

        vector<bool> undef(max_allele + 1, false);
        for (auto a : unique_alleles) {
                undef[a] = uniquekmers->is_undefined_allele(a);
        }

        const double one_third = 1.0 / 3.0;
        const double half = 0.5;

        for (size_t idx1 = 0; idx1 < unique_alleles.size(); ++idx1) {
                unsigned short a1 = unique_alleles[idx1];
                const bool a1_und = undef[a1];
                for (size_t idx2 = idx1; idx2 < unique_alleles.size(); ++idx2) {
                        unsigned short a2 = unique_alleles[idx2];
                        const bool a2_und = undef[a2];
                        double result = 1.0;
                        if (a1_und && a2_und) {
                                for (size_t i = 0; i < M; ++i) {
                                        result *= one_third * (P_table[i][0] + P_table[i][1] + P_table[i][2]);
                                }
                        } else if (a1_und || a2_und) {
                                for (size_t i = 0; i < M; ++i) {
                                        const int ec = (int)kmer_bit[i * row_stride + a1] + (int)kmer_bit[i * row_stride + a2];
                                        assert(ec < 2);
                                        result *= half * (P_table[i][ec] + P_table[i][ec + 1]);
                                }
                        } else {
                                for (size_t i = 0; i < M; ++i) {
                                        const int ec = (int)kmer_bit[i * row_stride + a1] + (int)kmer_bit[i * row_stride + a2];
                                        result *= P_table[i][ec];
                                }
                        }
                        this->state_to_prob[(size_t)a1 * stride + a2] = result;
                        this->state_to_prob[(size_t)a2 * stride + a1] = result;
                        if (result > 0) this->all_zeros = false;
                }
        }
}

double EmissionProbabilityComputer::get_emission_probability(unsigned short allele_id1, unsigned short allele_id2) const {
        if (this->all_zeros) return 1.0;
        return this->state_to_prob[(size_t)allele_id1 * ((size_t)this->max_allele + 1) + allele_id2];
}

double EmissionProbabilityComputer::compute_emission_probability(unsigned short allele_id1, unsigned short allele_id2, bool a1_undefined, bool a2_undefined){
        double result = 1.0;
        for (size_t i = 0; i < this->uniquekmers->size(); ++i){
                unsigned int expected_kmer_count = this->uniquekmers->kmer_on_allele(i, allele_id1) + this->uniquekmers->kmer_on_allele(i, allele_id2);
                if (a1_undefined && a2_undefined) {
                        result *= (1.0 / 3.0) * (this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(0) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(1) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(2));
                } else if (a1_undefined || a2_undefined) {
                        assert (expected_kmer_count < 2);
                        result *= 0.5 * (this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count) + this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count + 1));
                } else {
                        result *= this->probabilities->get_probability(this->uniquekmers->get_coverage(), this->uniquekmers->get_readcount_of(i)).get_probability_of(expected_kmer_count);
                }
        }
        return result;
}
