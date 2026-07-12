#ifndef B2SREADER_HPP
#define B2SREADER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <jellyfish/mer_dna.hpp>
#include "kmercounter.hpp"

/**
 * KmerCounter backed by a pre-computed back_to_sequences (b2s) count table.
 *
 * b2s (Peterlongo et al.) is run OUTSIDE PanGenie:
 *   back_to_sequences --in-kmers <panel_unique_kmers.fasta> \
 *                     --in-sequences <reads.fq> --out-kmers <table> -k 31
 * where the query fasta is the union of the `unique_kmers` (col4) and
 * `unique_kmers_overhang` (col5) columns of PanGenie-index's
 * `<prefix>_<chrom>_kmers.tsv.gz` files. Those are exactly the k-mers PanGenie
 * queries at genotype time (fill_read_kmercounts + compute_local_coverage), so
 * the table covers every getKmerAbundance() call.
 *
 * The table is `<kmer_string> <count>` (whitespace-separated, one per line);
 * b2s emits every query k-mer, canonical by default (a k-mer and its reverse
 * complement counted as one). This class reads that table into an in-memory
 * canonical map and answers getKmerAbundance / computeHistogram with semantics
 * byte-identical to KmcCounter (the abandoned stage6 KMC backend):
 *   * k-mer absent from reads (count 0) -> not stored, getKmerAbundance -> 0
 *   * k-mer present in reads            -> stored count (capped at 65535)
 * Count-0 rows are skipped on load (a missing key already returns 0), which
 * keeps the map to only the k-mers actually seen in the sample's reads.
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

class B2sReader : public KmerCounter {
public:
    /**
     * @param count_table path to the b2s `--out-kmers` table (`<kmer> <count>`
     *                    per line). May be gzip-compressed (name ending .gz).
     * @param kmer_size   k-mer length (must match PanGenie-index, e.g. 31).
     */
    B2sReader(const std::string& count_table, size_t kmer_size);
    ~B2sReader();

    size_t getKmerAbundance(std::string kmer) override;
    size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) override;
    size_t computeKmerCoverage(size_t genome_kmers) override;
    size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "") override;

private:
    std::unordered_map<jellyfish::mer_dna, uint16_t, MerDnaHash, MerDnaEqual> kmer_counts_;
    size_t kmer_size_;

    void load_table(const std::string& count_table);
};

#endif // B2SREADER_HPP
