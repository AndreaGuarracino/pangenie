#ifndef BIALLELICUNIQUEKMERS_HPP
#define BIALLELICUNIQUEKMERS_HPP

#include <array>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include "copynumber.hpp"
#include "kmerpath16.hpp"
#include "uniquekmers.hpp"
#include <cereal/types/polymorphic.hpp>
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>


/*
* Represents the set of unique kmers for a variant position.
*/


struct AlleleInfo16 {
	AlleleInfo16() : is_undefined(false), is_present(false) {}

	KmerPath16 kmer_path;
	bool is_undefined;
	// Presence used by the in-memory dense representation; not serialized.
	bool is_present;

	template <class Archive>
	void save(Archive& ar) const {
		ar(kmer_path, is_undefined);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(kmer_path, is_undefined);
		is_present = true;
	}
};


class BiallelicUniqueKmers : public UniqueKmers {
public:
	/**
	* @param variant_position genomic variant position
	* @param alleles defines which path (= index) covers each allele (= alleles[index])
	**/
	BiallelicUniqueKmers() = default;
	BiallelicUniqueKmers(size_t variant_position, std::vector<unsigned short>& alleles);

	size_t get_variant_position() const;
	void set_coverage(unsigned short local_coverage); 
	unsigned short get_coverage() const;

	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(unsigned short readcount, std::vector<unsigned short>& allele_ids);
	/** checks if kmer at index kmer_index is on path path_id **/
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	/** checks if kmer at index kmer_index is on allele allele_id 
	* Note: optimized for speed, no boundary checks performed.
	**/
	bool kmer_on_allele(size_t kmer_index, size_t allele_id) const;
	unsigned short get_readcount_of(size_t kmer_index);
	/** modify kmer count of an already inserted kmer **/
	void update_readcount(size_t kmer_index, unsigned short new_count);
	/** number of unique kmers **/
	size_t size() const;
	/** return number of paths **/
	unsigned short get_nr_paths() const;
	/** get all paths and alleles covering this position. If only_include, make sure to only output path_ids that are contained in only_include. **/
	void get_path_ids(std::vector<unsigned short>& paths, std::vector<unsigned short>& alleles, std::vector<unsigned short>* only_include = nullptr);
	/** get all unique alleles covered at this position **/
	void get_allele_ids(std::vector<unsigned short>& a);
	/** get only those unique alleles which are not undefined **/
	void get_defined_allele_ids(std::vector<unsigned short>& a);
	friend std::ostream& operator<< (std::ostream& stream, const BiallelicUniqueKmers& uk);
	/** returns a map which contains the number of unique kmers covering each allele **/
	std::map<unsigned short, int> kmers_on_alleles () const;
	/** returns the number of unique kmers on given allele */
	unsigned short kmers_on_allele(unsigned short allele_id) const;
	/** returns the number of read-supported kmers on given allele **/
	unsigned short present_kmers_on_allele(unsigned short allele_id) const;
	/** returns the fraction of read-supported kmers on given allele **/
	float fraction_present_kmers_on_allele(unsigned short allele_id) const;
	/** check whether allele is undefined **/
	bool is_undefined_allele (unsigned short allele_id) const;
	/** set allele to undefined **/
	void set_undefined_allele (unsigned short allele_id);
	/** look up allele covered by a path **/
	unsigned short get_allele(unsigned short path_id) const;
	/** update BiallelicUniqueKmers object by keeping only the paths provided **/
	void update_paths(std::vector<unsigned short>& path_ids);
	/** print kmer matrix (mainly for debugging) */
	void print_kmer_matrix(std::string chromosome) const;

	template<class Archive>
	void save(Archive& archive) const {
		std::map<bool, AlleleInfo16> serialized_alleles;
		for (size_t allele_id = 0; allele_id < alleles.size(); ++allele_id) {
			if (alleles[allele_id].is_present) {
				serialized_alleles.emplace(static_cast<bool>(allele_id), alleles[allele_id]);
			}
		}
		archive(variant_pos, local_coverage, current_index, kmer_to_count, serialized_alleles, path_to_allele);
	}

	template<class Archive>
	void load(Archive& archive) {
		std::map<bool, AlleleInfo16> serialized_alleles;
		archive(variant_pos, local_coverage, current_index, kmer_to_count, serialized_alleles, path_to_allele);
		alleles = {};
		for (const auto& entry : serialized_alleles) {
			alleles[entry.first] = entry.second;
			alleles[entry.first].is_present = true;
		}
	}

private:
	size_t variant_pos;
	float local_coverage;
	size_t current_index;
	std::vector<unsigned short> kmer_to_count;
	// Biallelic IDs are exactly 0 and 1; is_present distinguishes missing IDs.
	std::array<AlleleInfo16, 2> alleles;
	// defines which alleles are carried by each path (=index)
	std::vector<bool> path_to_allele;
//	friend class HaplotypeSampler;
	friend cereal::access;
};

#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/details/static_object.hpp>
CEREAL_REGISTER_TYPE(BiallelicUniqueKmers);
CEREAL_REGISTER_POLYMORPHIC_RELATION(UniqueKmers, BiallelicUniqueKmers)

# endif // BIALLELICUNIQUEKMERS_HPP
