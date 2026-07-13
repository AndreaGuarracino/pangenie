#include <algorithm>
#include <iostream>
#include <sstream>
#include "columnindexer.hpp"

using namespace std;

// NOTE: assumes that all UniqueKmers objects are covered by the same set of paths.
ColumnIndexer::ColumnIndexer (vector<shared_ptr<UniqueKmers>>* unique_kmers, vector<unsigned short>* only_paths)
	:unique_kmers(unique_kmers)
{
	size_t column_count = unique_kmers->size();
	this->variant_positions.reserve(column_count);
	for (size_t column_index = 0; column_index < column_count; ++ column_index) {
		vector<unsigned short> current_paths;
		vector<unsigned short> current_alleles;
		unique_kmers->at(column_index)->get_path_ids(current_paths, current_alleles, only_paths);
		unsigned short nr_paths = current_paths.size();	

		if (nr_paths == 0) {
			ostringstream oss;
			oss << "HMM::index_columns: column " << column_index << " is not covered by any paths.";
			throw runtime_error(oss.str());
		}
		if (column_index == 0) this->paths = current_paths;
		// check whether there are any non-reference alleles in panel
		bool all_absent = true;
		for (unsigned short i = 0; i < nr_paths; ++i) {
			if ((current_alleles[i] != 0) && (!unique_kmers->at(column_index)->is_undefined_allele(current_alleles[i])) ) all_absent = false;
		}
		if (!all_absent) {
			this->variant_positions.push_back(column_index);
			if (this->paths.size() != current_paths.size()) {
				throw runtime_error("ColumnIndexer::ColumnIndexer: varying number of paths across columns.");
			}
				if (current_paths == this->paths) {
					this->allele_by_column_path.insert(this->allele_by_column_path.end(), current_alleles.begin(), current_alleles.end());
				} else {
					// Preserve support for callers that provide the same paths in a different order.
					for (unsigned short path_id : this->paths) {
						auto path_it = find(current_paths.begin(), current_paths.end(), path_id);
						if (path_it == current_paths.end()) {
							throw runtime_error("ColumnIndexer::ColumnIndexer: path missing in column.");
						}
						this->allele_by_column_path.push_back(current_alleles[path_it - current_paths.begin()]);
					}
				}
		}
	}
}

pair<unsigned short,unsigned short> ColumnIndexer::get_path_ids_at (size_t position) const {
	if (position >= this->paths.size()*this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_path_ids_at: index out of bounds.");
	}
	size_t p_id1 = position % this->nr_paths();
	size_t p_id2 = (position / this->nr_paths()) % this->nr_paths();
	return pair<unsigned short,unsigned short> (p_id2, p_id1);
}
