#include "kmccounter.hpp"
#include "histogram.hpp"
#include "sequenceutils.hpp"
#include "kmerparser.hpp"
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <zlib.h>

using namespace std;

static string make_tmp_dir() {
    string tmpl = "/tmp/pangenie_kmc_XXXXXX";
    vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (mkdtemp(buf.data()) == nullptr) {
        throw runtime_error("KmcCounter: mkdtemp failed");
    }
    return string(buf.data());
}

static void remove_dir(const string& path) {
    string cmd = "rm -rf '" + path + "'";
    int rc = system(cmd.c_str());
    (void)rc;
}

KmcCounter::KmcCounter(const string& readfile,
                       size_t kmer_size,
                       size_t nr_threads,
                       const vector<string>& panel_kmer_tsv_files,
                       uint64_t /*hash_size*/)
    : kmer_size_(kmer_size)
{
    jellyfish::mer_dna::k(kmer_size);
    tmp_dir_ = make_tmp_dir();
    build_panel_set_and_run(readfile, kmer_size, nr_threads, panel_kmer_tsv_files);
}

KmcCounter::~KmcCounter() {
    if (!tmp_dir_.empty()) remove_dir(tmp_dir_);
}

void KmcCounter::build_panel_set_and_run(const string& readfile,
                                         size_t kmer_size,
                                         size_t nr_threads,
                                         const vector<string>& panel_kmer_tsv_files)
{
    // STAGE 1: walk the per-chromosome _kmers.tsv.gz files to build the canonical
    // panel kmer set. These are the only kmers PanGenie will ever query (both the
    // unique-kmer-per-variant set used in fill_read_kmercounts and the flanking
    // kmers used in compute_local_coverage).
    unordered_set<jellyfish::mer_dna, MerDnaHash, MerDnaEqual> panel_set;
    panel_set.reserve(64ULL * 1024 * 1024);  // initial guess; HPRC v1 panel is ~16M-100M kmers

    cerr << "[KmcCounter] building panel kmer set from "
         << panel_kmer_tsv_files.size() << " kmer-tsv file(s)" << endl;

    for (const auto& path : panel_kmer_tsv_files) {
        gzFile fh = gzopen(path.c_str(), "rb");
        if (!fh) throw runtime_error("KmcCounter: cannot open panel kmer file: " + path);
        const int buffer_size = 1024;
        char buf[buffer_size];
        string line;
        while (gzgets(fh, buf, buffer_size) != nullptr) {
            line += buf;
            if (line.back() != '\n') continue;
            line.pop_back();
            string chrom;
            size_t start;
            vector<string> kmers;
            vector<string> flanking_kmers;
            bool is_header = false;
            parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);
            line.clear();
            if (is_header) continue;
            for (auto& k : kmers) {
                if (k.size() != kmer_size) continue;
                jellyfish::mer_dna mer(k);
                mer.canonicalize();
                panel_set.insert(std::move(mer));
            }
            for (auto& k : flanking_kmers) {
                if (k.size() != kmer_size) continue;
                jellyfish::mer_dna mer(k);
                mer.canonicalize();
                panel_set.insert(std::move(mer));
            }
        }
        gzclose(fh);
    }
    cerr << "[KmcCounter] panel kmer set has " << panel_set.size() << " unique canonical kmers" << endl;

    // STAGE 2: run KMC on the read inputs. KMC accepts .gz files directly and a
    // whitespace-separated list (we rewrite to an @file-list when more than one).
    string kmc_input;
    if (!readfile.empty() && readfile[0] == '@') {
        kmc_input = readfile;
    } else {
        istringstream iss(readfile);
        vector<string> paths;
        string tok;
        while (iss >> tok) paths.push_back(tok);
        if (paths.size() <= 1) {
            kmc_input = "'" + readfile + "'";
        } else {
            string list_path = tmp_dir_ + "/inputs.lst";
            ofstream out(list_path);
            for (auto& p : paths) out << p << "\n";
            out.close();
            kmc_input = "@" + list_path;
        }
    }

    string db_prefix = tmp_dir_ + "/db";
    ostringstream cmd;
    cmd << "kmc -k" << kmer_size
        << " -fq -ci1 -cs65535"
        << " -t" << nr_threads
        << " -m24"
        << " " << kmc_input
        << " '" << db_prefix << "'"
        << " '" << tmp_dir_ << "' >&2";
    cerr << "[KmcCounter] running: " << cmd.str() << endl;
    int rc = system(cmd.str().c_str());
    if (rc != 0) {
        throw runtime_error("KmcCounter: kmc binary failed (exit " + to_string(rc) + ")");
    }

    // STAGE 3: walk KMC output and insert ONLY kmers present in panel_set.
    CKMCFile kmc_db;
    if (!kmc_db.OpenForListing(db_prefix)) {
        throw runtime_error("KmcCounter: CKMCFile::OpenForListing failed for " + db_prefix);
    }
    CKMCFileInfo info;
    kmc_db.Info(info);
    cerr << "[KmcCounter] KMC counted " << info.total_kmers
         << " distinct read kmers (k=" << info.kmer_length << "); filtering to panel..." << endl;

    CKmerAPI kmer((uint32)kmer_size);
    string kmer_str(kmer_size, 'N');
    uint64 count = 0;
    kmer_counts_.reserve(panel_set.size());
    size_t kept = 0;
    while (kmc_db.ReadNextKmer(kmer, count)) {
        kmer.to_string((char*)kmer_str.data());
        jellyfish::mer_dna mer(kmer_str);
        mer.canonicalize();
        if (panel_set.find(mer) == panel_set.end()) continue;
        uint16_t capped = (count > 65535) ? (uint16_t)65535 : (uint16_t)count;
        kmer_counts_.emplace(std::move(mer), capped);
        ++kept;
    }
    kmc_db.Close();
    cerr << "[KmcCounter] kept " << kept << " of " << panel_set.size()
         << " panel kmers in the in-memory hash ("
         << (panel_set.size() ? (100.0 * kept / panel_set.size()) : 0.0) << "%)" << endl;
}

size_t KmcCounter::getKmerAbundance(string kmer) {
    jellyfish::mer_dna mer(kmer);
    mer.canonicalize();
    auto it = kmer_counts_.find(mer);
    return (it == kmer_counts_.end()) ? 0 : it->second;
}

size_t KmcCounter::getKmerAbundance(jellyfish::mer_dna jelly_kmer) {
    jelly_kmer.canonicalize();
    auto it = kmer_counts_.find(jelly_kmer);
    return (it == kmer_counts_.end()) ? 0 : it->second;
}

size_t KmcCounter::computeKmerCoverage(size_t genome_kmers) {
    double result = 0.0;
    double genome = 1.0 * genome_kmers;
    for (const auto& kv : kmer_counts_) {
        result += (1.0 * kv.second) / genome;
    }
    return (size_t) ceil(result);
}

size_t KmcCounter::computeHistogram(size_t max_count, bool largest_peak, string filename) {
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
            ss << "KmcCounter::computeHistogram: File " << filename << " cannot be created." << endl;
            throw runtime_error(ss.str());
        }
        histofile << "parameters\t" << kmer_coverage_estimate/2.0 << '\t' << kmer_coverage_estimate << endl;
        histofile.close();
    }
    return kmer_coverage_estimate;
}
