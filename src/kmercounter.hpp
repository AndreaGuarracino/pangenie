#ifndef KMERCOUNTER_HPP
#define KMERCOUNTER_HPP

#include <map>
#include <vector>
#include <string>
#include <string_view>
#include <span>
#include <stdexcept>
#include <jellyfish/mer_dna.hpp>

class KmerCounter {
public:
	/** get the abundance of given kmer (string) **/
	virtual size_t getKmerAbundance(std::string_view kmer) = 0;
	virtual void getKmerAbundances(std::span<const std::string_view> kmers, std::span<size_t> counts) {
		if (kmers.size() != counts.size()) throw std::invalid_argument("KmerCounter::getKmerAbundances: size mismatch.");
		for (size_t i = 0; i < kmers.size(); ++i) counts[i] = getKmerAbundance(kmers[i]);
	}

	/** get the abundance of given kmer (jellyfish kmer) **/
	virtual size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) = 0;

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	virtual size_t computeKmerCoverage(size_t genome_kmers) = 0;

	/** computes kmer abundance histogram and returns the three highest peaks **/
	virtual size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "") = 0;

	virtual ~KmerCounter() {} ;
};
#endif // KMERCOUNTER_HPP
