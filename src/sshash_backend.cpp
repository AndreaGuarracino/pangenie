#include "sshash_backend.hpp"

// SSHash's own "include the .cpp" pattern (mirrors tools/sshash.cpp). These carry the
// template method definitions; instantiated for dictionary_type when we use it below.
// Requires -I <sshash_root> and the pthash/bits/essentials/fastmod/cityhash include dirs.
#include "include/dictionary_types.hpp"
#include "src/dictionary.cpp"
#include "src/query.cpp"
#include "src/info.cpp"
// cityhash (the minimizer hasher) is a plain .cpp — compile it into this TU.
#include "external/cityhash/cityhash.cpp"

#include <essentials.hpp>  // essentials::load

// Edit E: streaming read counter lives in THIS TU so streaming_query.lookup and
// the count bump inline together in the hot loop (no cross-TU call per k-mer).
#include "include/streaming_query.hpp"
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <atomic>
#include <cstdint>
#include <iostream>

namespace {
struct Impl {
	sshash::dictionary_type dict;
};

typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sd_parser_t;

// Aggregated streaming stats (across threads) for one count_reads_streaming call.
std::atomic<uint64_t> g_num_searches{0};    // full minimizer+MPHF probes
std::atomic<uint64_t> g_num_extensions{0};  // cheap in-unitig extensions (the win)

// Per-thread streaming walker: pulls parser jobs, resets the streaming_query at
// each job boundary (matches jellyfish's filled_=0 reset), and walks every
// k-window in the job buffer. streaming_query.lookup handles validity ('N'
// between reads -> reset) and extension internally; its returned kmer_id is
// identical to dictionary::lookup, so counts are byte-identical.
template <bool CANON>
class StreamCountScanner : public jellyfish::thread_exec {
	jellyfish::stream_manager<char**> streams_;
	sd_parser_t parser_;
	const sshash::dictionary_type* dict_;
	uint32_t* counts_;
	const uint64_t k_;
public:
	StreamCountScanner(int nb_threads, char** fb, char** fe,
	                   const sshash::dictionary_type* dict, uint32_t* counts)
		: streams_(fb, fe)
		, parser_((uint16_t)dict->k(), streams_.nb_streams(), 3 * nb_threads, 64 * 1024, streams_)
		, dict_(dict), counts_(counts), k_(dict->k())
	{}

	virtual void start(int /*thid*/) {
		sshash::streaming_query<sshash::dictionary_type, CANON> q(dict_);
		typename sd_parser_t::job job(parser_);
		for (; !job.is_empty(); job.next()) {
			const char* s = job->start;
			const char* e = job->end;
			q.reset();
			for (const char* p = s; p + (long)k_ <= e; ++p) {
				auto r = q.lookup(p);
				if (r.kmer_id != sshash::constants::invalid_uint64) {
					std::atomic_ref<uint32_t> ref(counts_[r.kmer_id]);
					ref.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		// num_searches()/num_extensions() are cumulative over this query's lifetime
		// (reset() does NOT clear them), so aggregate once per thread, after all jobs.
		g_num_searches.fetch_add(q.num_searches(), std::memory_order_relaxed);
		g_num_extensions.fetch_add(q.num_extensions(), std::memory_order_relaxed);
	}
};
}  // namespace

SSHashDict::SSHashDict() : impl_(new Impl()) {}
SSHashDict::~SSHashDict() { delete static_cast<Impl*>(impl_); }

bool SSHashDict::load(const std::string& path) {
	auto* d = &static_cast<Impl*>(impl_)->dict;
	essentials::load(*d, path.c_str());
	return d->num_kmers() > 0;
}

uint64_t SSHashDict::num_kmers() const {
	return static_cast<Impl*>(impl_)->dict.num_kmers();
}

uint64_t SSHashDict::lookup(const char* kmer) const {
	auto r = static_cast<const Impl*>(impl_)->dict.lookup(kmer, /*check_reverse_complement=*/true);
	return (r.kmer_id == sshash::constants::invalid_uint64) ? NOT_FOUND : r.kmer_id;
}

void SSHashDict::count_reads_streaming(char** file_begin, char** file_end,
                                       int nr_threads, uint32_t* counts) const {
	const sshash::dictionary_type* dict = &static_cast<const Impl*>(impl_)->dict;
	g_num_searches.store(0, std::memory_order_relaxed);
	g_num_extensions.store(0, std::memory_order_relaxed);
	if (dict->canonical()) {
		StreamCountScanner<true> sc(nr_threads, file_begin, file_end, dict, counts);
		sc.exec_join(nr_threads);
	} else {
		StreamCountScanner<false> sc(nr_threads, file_begin, file_end, dict, counts);
		sc.exec_join(nr_threads);
	}
	uint64_t se = g_num_searches.load(std::memory_order_relaxed);
	uint64_t ex = g_num_extensions.load(std::memory_order_relaxed);
	double tot = (double)(se + ex);
	std::cerr << "[count-phase] stream_searches=" << se << " stream_extensions=" << ex
	          << " extension_frac=" << (tot > 0 ? ex / tot : 0.0) << '\n';
}
