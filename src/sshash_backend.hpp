#ifndef SSHASH_BACKEND_HPP
#define SSHASH_BACKEND_HPP

#include <string>
#include <cstdint>

/**
 * Edit D1: thin pimpl wrapper around an SSHash dictionary.
 *
 * Hides all of SSHash's templates behind a plain interface so only ONE translation
 * unit (sshash_backend.cpp) compiles the heavy SSHash headers, and the rest of
 * PanGenie stays clean. The dictionary is built OFFLINE (ggcat -> sshash build) and
 * loaded here; it maps each canonical panel k-mer to a dense id in [0, num_kmers).
 */
class SSHashDict {
public:
	static constexpr uint64_t NOT_FOUND = ~0ULL;

	SSHashDict();
	~SSHashDict();
	SSHashDict(const SSHashDict&) = delete;
	SSHashDict& operator=(const SSHashDict&) = delete;

	// load a serialized .sshash dictionary; returns true on success
	bool load(const std::string& path);

	uint64_t num_kmers() const;

	// look up a k-length k-mer string (canonicalized internally); returns the
	// dense kmer-id, or NOT_FOUND if the k-mer is not in the panel.
	uint64_t lookup(const char* kmer) const;

	// Edit E: streaming read counter. Parses the read file(s) with jellyfish's
	// mer_overlap_sequence_parser (identical file/read/seam/N handling as the
	// per-kmer path) and, per maximal valid run, uses an SSHash streaming_query
	// to look up each consecutive read k-mer -- extending within a unitig (a
	// base compare) instead of a fresh minimizer+MPHF probe per k-mer. Bumps
	// counts[kmer_id] atomically. Byte-identical to per-kmer lookup by
	// construction: streaming_query.lookup returns the same kmer_id as lookup().
	void count_reads_streaming(char** file_begin, char** file_end,
	                           int nr_threads, uint32_t* counts) const;

private:
	void* impl_;
};

#endif // SSHASH_BACKEND_HPP
