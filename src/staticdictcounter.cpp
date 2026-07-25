#include "staticdictcounter.hpp"
#include "histogram.hpp"
#include "sequenceutils.hpp"      // compute_kmer_coverage
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <stdexcept>

using namespace std;

namespace {
vector<char*> sd_to_args(const string& readfile) {
	vector<char*> args;
	istringstream iss(readfile);
	string token;
	while (iss >> token) {
		char* arg = new char[token.size() + 1];
		copy(token.begin(), token.end(), arg);
		arg[token.size()] = '\0';
		args.push_back(arg);
	}
	args.push_back(0);
	return args;
}
}  // namespace

StaticDictCounter::StaticDictCounter(string readfile, string sshash_path,
                                     size_t kmer_size, size_t nr_threads, uint64_t /*hash*/)
	: nr_threads_(nr_threads), k_(kmer_size)
{
	jellyfish::mer_dna::k(kmer_size);
	if (sshash_path.empty())
		throw runtime_error("StaticDictCounter requires an SSHash dictionary path (build it with PanGenie-index -S)");
	auto t0 = chrono::steady_clock::now();
	if (!sshash_.load(sshash_path))
		throw runtime_error(string("StaticDictCounter: failed to load SSHash dictionary '") + sshash_path + "'");
	counts_.assign(sshash_.num_kmers(), 0u);
	double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
	cerr << "[count-phase] sshash_load_s=" << dt << " panel_keys=" << sshash_.num_kmers()
	     << "  (per-sample PRIME eliminated)\n";
	count_reads(readfile);
}

void StaticDictCounter::count_reads(const string& readfile) {
	auto t0 = chrono::steady_clock::now();
	vector<char*> args = sd_to_args(readfile);
	// Streaming SSHash query: extends within unitigs instead of a fresh
	// minimizer+MPHF probe per read k-mer. Byte-identical to a per-kmer lookup.
	sshash_.count_reads_streaming(&args[0], (&args[0]) + 1, (int)nr_threads_, counts_.data());
	for (auto* a : args) delete[] a;
	double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
	cerr << "[count-phase] dict_count_reads_s=" << dt << '\n';
}

size_t StaticDictCounter::getKmerAbundance(string_view kmer) {
	const size_t k = jellyfish::mer_dna::k();
	if (kmer.size() < k) return 0;
	char buf[64];
	memcpy(buf, kmer.data(), k);
	buf[k] = '\0';
	uint64_t id = sshash_.lookup(buf);
	return (id == SSHashDict::NOT_FOUND) ? 0 : (size_t)counts_[id];
}

size_t StaticDictCounter::getKmerAbundance(jellyfish::mer_dna jelly_kmer) {
	char buf[64];
	jelly_kmer.to_str(buf);
	buf[jellyfish::mer_dna::k()] = '\0';
	uint64_t id = sshash_.lookup(buf);
	return (id == SSHashDict::NOT_FOUND) ? 0 : (size_t)counts_[id];
}

size_t StaticDictCounter::computeKmerCoverage(size_t genome_kmers) {
	double result = 0.0;
	double genome = 1.0 * genome_kmers;
	for (uint32_t c : counts_) result += (1.0 * c) / genome;
	return (size_t)ceil(result);
}

// Mirrors JellyfishCounter::computeHistogram exactly, scanning the dense count array.
size_t StaticDictCounter::computeHistogram(size_t max_count, bool largest_peak, string filename) {
	Histogram histogram(max_count);
	const size_t worker_count = max<size_t>(1, nr_threads_);
	vector<vector<size_t>> partial_histograms(worker_count, vector<size_t>(max_count + 1, 0));
	vector<thread> workers;
	workers.reserve(worker_count);
	const size_t n = counts_.size();
	auto t0 = chrono::steady_clock::now();
	for (size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
		workers.emplace_back([this, worker_id, worker_count, max_count, n, &partial_histograms]() {
			auto& bins = partial_histograms[worker_id];
			const size_t lo = n / worker_count * worker_id + std::min(n % worker_count, worker_id);
			const size_t hi = n / worker_count * (worker_id + 1) + std::min(n % worker_count, worker_id + 1);
			for (size_t i = lo; i < hi; ++i) {
				size_t value = counts_[i];
				if ((value > 0) && (value <= max_count)) ++bins[value];
			}
		});
	}
	for (auto& worker : workers) worker.join();
	double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
	size_t nonzero = 0;
	for (size_t value = 1; value <= max_count; ++value) {
		size_t count = 0;
		for (const auto& bins : partial_histograms) count += bins[value];
		histogram.add_count(value, count);
		nonzero += count;
	}
	cerr << "[count-phase] dict_histogram_scan_s=" << dt
	     << " panel_kmers_1.._max=" << nonzero
	     << " entries_scanned=" << n << '\n';

	if (filename != "") histogram.write_to_file(filename);
	histogram.smooth_histogram();
	vector<size_t> peak_ids, peak_values;
	histogram.find_peaks(peak_ids, peak_values);
	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, largest_peak);
	if (filename != "") {
		ofstream histofile;
		histofile.open(filename, ios::app);
		if (!histofile.good()) {
			stringstream ss;
			ss << "StaticDictCounter::computeHistogram: File " << filename << " cannot be created." << endl;
			throw runtime_error(ss.str());
		}
		histofile << "parameters\t" << kmer_coverage_estimate / 2.0 << '\t' << kmer_coverage_estimate << '\n';
		histofile.close();
	}
	return kmer_coverage_estimate;
}
