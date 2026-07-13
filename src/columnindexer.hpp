#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>
#include <memory>
#include <span>
#include <stdexcept>
#include "uniquekmers.hpp"

/** 
* Keep track of paths and alleles in a column
**/

class ColumnIndexer {
public:
	ColumnIndexer(std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers, std::vector<unsigned short>* only_paths);
	/** get variant ID corresponding to column index **/
	size_t get_variant_id(size_t column_index) const {
		if (column_index >= variant_positions.size()) throw std::runtime_error("ColumnIndexer::get_variant_id: column index does not exist.");
		return variant_positions[column_index];
	}
	/** get the number of columns **/
	size_t size() const noexcept { return variant_positions.size(); }
	/** number of paths **/
	unsigned short nr_paths() const noexcept { return paths.size(); }
	/** get path_id at a given index **/
	unsigned short get_path(unsigned short path_index) const {
		if (path_index >= paths.size()) throw std::runtime_error("ColumnIndexer::get_path: path_index does not exist.");
		return paths[path_index];
	}
	/** get the allele covered by a path in a column **/
	unsigned short get_allele(unsigned short path_index, size_t column_index) const {
		if (path_index >= paths.size()) throw std::runtime_error("ColumnIndexer::get_allele: path_index does not exist.");
		if (column_index >= variant_positions.size()) throw std::runtime_error("ColumnIndex::get_allele: column_index does not exist.");
		return allele_by_column_path[column_index * paths.size() + path_index];
	}
	std::span<const unsigned short> get_alleles(size_t column_index) const {
		if (column_index >= variant_positions.size()) throw std::runtime_error("ColumnIndex::get_alleles: column_index does not exist.");
		return {allele_by_column_path.data() + column_index * paths.size(), paths.size()};
	}
	std::span<const unsigned short> get_paths() const noexcept { return paths; }
	/** get path ids corresponding to an index inside of a column **/
	std::pair<unsigned short,unsigned short> get_path_ids_at (size_t position) const;

private:
	std::vector<size_t> variant_positions;
	std::vector<unsigned short> paths;
	std::vector<unsigned short> allele_by_column_path;
	std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers;
};

#endif // COLUMNINDEXER_HPP
