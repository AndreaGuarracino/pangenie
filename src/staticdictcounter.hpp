#ifndef STATICDICTCOUNTER_HPP
#define STATICDICTCOUNTER_HPP

#include <vector>
#include <cstdint>
#include <string>
#include <string_view>
#include <jellyfish/mer_dna.hpp>
#include "kmercounter.hpp"
#include "sshash_backend.hpp"

/**
 * SSHash-backed static-dictionary read counter.
 *
 * Replaces the Jellyfish PRIME+UPDATE design. An immutable SSHash dictionary
 * (canonical panel k-mer -> dense id in [0, num_kmers)) is built ONCE per panel
 * by PanGenie-index -S (ggcat -> sshash) and loaded here from <index_prefix>.sshash,
 * so there is NO per-sample PRIME. Read k-mers are streamed through the dictionary
 * and counted into a dense uint32 array. Histogram and coverage are computed over
 * that array, mirroring JellyfishCounter exactly (byte-identical genotypes).
 */
class StaticDictCounter : public KmerCounter {
public:
	// sshash_path: path to the prebuilt .sshash dictionary (PanGenie-index -S writes
	// <prefix>.sshash). Empty path is an error.
	StaticDictCounter(std::string readfile, std::string sshash_path,
	                  size_t kmer_size, size_t nr_threads, uint64_t hash = 0);
	~StaticDictCounter() override = default;

	size_t getKmerAbundance(std::string_view kmer) override;
	size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) override;
	size_t computeKmerCoverage(size_t genome_kmers) override;
	size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "") override;

private:
	size_t nr_threads_;
	size_t k_;
	SSHashDict sshash_;              // immutable panel dictionary (loaded from file)
	std::vector<uint32_t> counts_;   // dense read counts, size = sshash_.num_kmers()

	void count_reads(const std::string& readfile);
};

#endif // STATICDICTCOUNTER_HPP
