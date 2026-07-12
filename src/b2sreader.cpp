#include "b2sreader.hpp"
#include "histogram.hpp"
#include "sequenceutils.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <math.h>
#include <zlib.h>

using namespace std;

B2sReader::B2sReader(const string& count_table, size_t kmer_size)
    : kmer_size_(kmer_size)
{
    // mer_dna is a global-state singleton: set k before constructing any mer_dna.
    jellyfish::mer_dna::k(kmer_size);
    load_table(count_table);
}

B2sReader::~B2sReader() {}

void B2sReader::load_table(const string& count_table) {
    // gzopen transparently reads both plain-text and gzip files.
    gzFile fh = gzopen(count_table.c_str(), "rb");
    if (!fh) throw runtime_error("B2sReader: cannot open b2s count table: " + count_table);

    cerr << "[B2sReader] loading b2s count table " << count_table << endl;
    kmer_counts_.reserve(8ULL * 1024 * 1024);

    const int buffer_size = 4096;
    char buf[buffer_size];
    string line;
    size_t rows = 0, kept = 0, skipped_zero = 0, skipped_len = 0;

    while (gzgets(fh, buf, buffer_size) != nullptr) {
        line += buf;
        if (line.empty() || line.back() != '\n') continue;  // wait for full line
        line.pop_back();
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') { line.clear(); continue; }

        // b2s --out-kmers format: "<kmer_string> <count>" (whitespace-separated).
        istringstream iss(line);
        string kmer;
        long long count = 0;
        iss >> kmer >> count;
        line.clear();
        rows++;

        if (kmer.size() != kmer_size_) { skipped_len++; continue; }
        if (count <= 0) { skipped_zero++; continue; }  // missing key already returns 0

        jellyfish::mer_dna mer(kmer);
        mer.canonicalize();
        uint16_t capped = (count > 65535) ? (uint16_t)65535 : (uint16_t)count;
        kmer_counts_.emplace(std::move(mer), capped);
        kept++;
    }
    gzclose(fh);

    cerr << "[B2sReader] read " << rows << " rows; kept " << kept
         << " nonzero canonical kmers (skipped " << skipped_zero << " zero-count, "
         << skipped_len << " wrong-length)" << endl;
    if (kept == 0) {
        throw runtime_error("B2sReader: no nonzero kmers loaded from " + count_table
                            + " (wrong format, wrong k, or empty overlap?)");
    }
}

size_t B2sReader::getKmerAbundance(string kmer) {
    jellyfish::mer_dna mer(kmer);
    mer.canonicalize();
    auto it = kmer_counts_.find(mer);
    return (it == kmer_counts_.end()) ? 0 : it->second;
}

size_t B2sReader::getKmerAbundance(jellyfish::mer_dna jelly_kmer) {
    jelly_kmer.canonicalize();
    auto it = kmer_counts_.find(jelly_kmer);
    return (it == kmer_counts_.end()) ? 0 : it->second;
}

size_t B2sReader::computeKmerCoverage(size_t genome_kmers) {
    double result = 0.0;
    double genome = 1.0 * genome_kmers;
    for (const auto& kv : kmer_counts_) {
        result += (1.0 * kv.second) / genome;
    }
    return (size_t) ceil(result);
}

size_t B2sReader::computeHistogram(size_t max_count, bool largest_peak, string filename) {
    Histogram histogram(max_count);
    for (const auto& kv : kmer_counts_) {
        if (kv.second > 0) histogram.add_value(kv.second);
    }
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
            ss << "B2sReader::computeHistogram: File " << filename << " cannot be created." << endl;
            throw runtime_error(ss.str());
        }
        histofile << "parameters\t" << kmer_coverage_estimate/2.0 << '\t' << kmer_coverage_estimate << endl;
        histofile.close();
    }
    return kmer_coverage_estimate;
}
