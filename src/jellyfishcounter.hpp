#ifndef JELLYFISHCOUNTER_HPP
#define JELLYFISHCOUNTER_HPP

#include <map>
#include <vector>
#include <string>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include "kmercounter.hpp"
#include "sequenceutils.hpp"

/**
* Counts Kmers in DNA-sequences (given in FASTQ-format) using jellyfish.
 **/


enum OPERATION { COUNT, PRIME, UPDATE };

// using example from: Jellyfish-2/examples/jf_count_dump/jf_count_dump.cc
typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>	mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>>	sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna>	mer_iterator_type;

class mer_counter : public jellyfish::thread_exec {
	mer_hash_type& mer_hash_;
	jellyfish::stream_manager<char**> streams_;
	sequence_parser_type parser_;
  const bool canonical_;
	OPERATION op_;

public:
	mer_counter(int nb_threads, mer_hash_type& mer_hash,
	char** file_begin, char** file_end,
	bool canonical, OPERATION op)
	: mer_hash_(mer_hash)
	, streams_(file_begin, file_end)
	, parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 64 * 1024, streams_)
	, canonical_(canonical)
	, op_(op)
{ }

	virtual void start(int /*thid*/) {
		mer_iterator_type mers(parser_, canonical_);

		switch(op_) {
			case COUNT:
				for( ; mers; ++mers)
					mer_hash_.add(*mers, 1);
				break;

			case PRIME:
				for( ; mers; ++mers)
					mer_hash_.set(*mers);
				break;

			case UPDATE:
				jellyfish::mer_dna tmp;
				for( ; mers; ++mers)
					mer_hash_.update_add(*mers, 1, tmp);
				break;
		}

		mer_hash_.done();
	}
};


// ---- Edit B: pluggable read scanner ----------------------------------------
// Reuses jellyfish's EXACT parser + mer_iterator (canonicalization, N-handling,
// read boundaries identical to mer_counter above), so the canonical k-mer stream
// is byte-identical by construction. Only the per-k-mer SINK changes. This is the
// single hook the membership filter (Edit C) and static dictionary (Edit D) plug into.
// Sink must provide: operator()(int thid, const jellyfish::mer_dna&) and finish(int thid).
template<class Sink>
class mer_scanner : public jellyfish::thread_exec {
	jellyfish::stream_manager<char**> streams_;
	sequence_parser_type parser_;
	const bool canonical_;
	Sink& sink_;
public:
	mer_scanner(int nb_threads, char** file_begin, char** file_end, bool canonical, Sink& sink)
	: streams_(file_begin, file_end)
	, parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 64 * 1024, streams_)
	, canonical_(canonical)
	, sink_(sink)
	{ }

	virtual void start(int thid) {
		mer_iterator_type mers(parser_, canonical_);
		for( ; mers; ++mers) sink_(thid, *mers);
		sink_.finish(thid);
	}
};

// Validation sink for Edit B: reproduces mer_counter's UPDATE exactly
// (per-thread scratch mer_dna; per-thread hash done()).
struct HashUpdateSink {
	mer_hash_type& hash_;
	std::vector<jellyfish::mer_dna> tmp_;
	HashUpdateSink(mer_hash_type& hash, int nb_threads) : hash_(hash), tmp_(nb_threads) {}
	inline void operator()(int thid, const jellyfish::mer_dna& m) { hash_.update_add(m, 1, tmp_[thid]); }
	inline void finish(int /*thid*/) { hash_.done(); }
};
// ---------------------------------------------------------------------------


class JellyfishCounter : public KmerCounter {
public:
	/** 
	* @param readfile name of the FASTQ-files containing reads
	* @param *params parameters for GATB-Kmercounter
	* @param name of the output file
	**/
	JellyfishCounter(std::string readfile, size_t kmer_size, size_t nr_threads = 1, uint64_t hash = 3000000000);

	/** 
	* @param readfile name of the FASTQ-files containing reads
	* @param kmerfile only count kmers contained in sequences given in these FASTQ-files
	* @param *params parameters for GATB-Kmercounter
	* @param name of the output file
	**/
	JellyfishCounter (std::string readfile, std::vector<std::string> kmerfiles, size_t kmer_size, size_t nr_threads = 1, uint64_t hash = 3000000000);

	~JellyfishCounter();
	
	/** get the abundance of given kmer (string) **/
	size_t getKmerAbundance(std::string_view kmer) override;
	void getKmerAbundances(std::span<const std::string_view> kmers, std::span<size_t> counts) override;

	/** get the abundance of given kmer (jellyfish kmer) **/
	size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) override;

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	size_t computeKmerCoverage(size_t genome_kmers);

	/** computes kmer abundance histogram and returns the three highest peaks **/
	size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "");

private:
	mer_hash_type* jellyfish_hash;
	size_t nr_threads;
};
#endif // JELLYFISHCOUNTER_HPP
