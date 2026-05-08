#ifndef KMCCOUNTER_HPP
#define KMCCOUNTER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <jellyfish/mer_dna.hpp>
#include "kmercounter.hpp"

/**
 * KmerCounter backed by KMC3, with panel-kmer filtering.
 *
 * Counts k-mers in `readfile` using KMC, then filters the KMC output against
 * the union of all kmers + flanking_kmers parsed from `panel_kmer_tsv_files`.
 * Only kmers present in the panel set are kept in the in-memory hash; the rest
 * of the read kmer space (typically billions of singleton sequencing errors)
 * is discarded.
 *
 * This mirrors JellyfishCounter's PRIME+UPDATE semantics for count_only_graph=true:
 *   * panel kmers absent from reads -> getKmerAbundance returns 0
 *   * panel kmers present in reads  -> getKmerAbundance returns the read count
 *   * non-panel kmers               -> never queried, so storage discarded
 *
 * KMC reads .gz files natively, so on AoU we can pass R1.fq.gz + R2.fq.gz
 * directly via whitespace-separated paths in `readfile` (no pigz step).
 */
struct MerDnaHash {
    size_t operator()(const jellyfish::mer_dna& m) const noexcept {
        // For k <= 32, nb_words() == 1; data()[0] uniquely identifies the kmer.
        return std::hash<uint64_t>{}(m.data()[0]);
    }
};

struct MerDnaEqual {
    bool operator()(const jellyfish::mer_dna& a, const jellyfish::mer_dna& b) const noexcept {
        return a == b;
    }
};

class KmcCounter : public KmerCounter {
public:
    /**
     * @param readfile path (or whitespace-separated paths or @file-list) for KMC input.
     *                 KMC accepts .gz files natively.
     * @param kmer_size k-mer length (e.g. 31).
     * @param nr_threads threads for KMC.
     * @param panel_kmer_tsv_files list of `<prefix>_<chrom>_kmers.tsv.gz` files written
     *                 by PanGenie-index. Only kmers present in these files are kept
     *                 in the in-memory hash.
     */
    KmcCounter(const std::string& readfile,
               size_t kmer_size,
               size_t nr_threads,
               const std::vector<std::string>& panel_kmer_tsv_files,
               uint64_t hash_size = 3000000000);
    ~KmcCounter();

    size_t getKmerAbundance(std::string kmer) override;
    size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) override;
    size_t computeKmerCoverage(size_t genome_kmers) override;
    size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "") override;

private:
    std::unordered_map<jellyfish::mer_dna, uint16_t, MerDnaHash, MerDnaEqual> kmer_counts_;
    std::string tmp_dir_;
    size_t kmer_size_;

    void build_panel_set_and_run(const std::string& readfile,
                                 size_t kmer_size,
                                 size_t nr_threads,
                                 const std::vector<std::string>& panel_kmer_tsv_files);
};

#endif
