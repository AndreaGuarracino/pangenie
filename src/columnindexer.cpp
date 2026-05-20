#include <iostream>
#include <sstream>
#include "columnindexer.hpp"

using namespace std;

// NOTE: assumes that all UniqueKmers objects are covered by the same set of paths.
ColumnIndexer::ColumnIndexer (vector<shared_ptr<UniqueKmers>>* unique_kmers, vector<unsigned short>* only_paths)
	:unique_kmers(unique_kmers)
{
	size_t column_count = unique_kmers->size();
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
			size_t nr_paths_size_t = this->paths.size();
			for (size_t path_index = 0; path_index < nr_paths_size_t; ++path_index) {
				unsigned short path_id = this->paths[path_index];
				bool found = false;
				for (size_t i = 0; i < current_paths.size(); ++i) {
					if (current_paths[i] == path_id) {
						this->allele_by_column_path.push_back(current_alleles[i]);
						found = true;
						break;
					}
				}
				if (!found) {
					throw runtime_error("ColumnIndexer::ColumnIndexer: path missing in column.");
				}
			}
		}
	}
}

size_t ColumnIndexer::get_variant_id(size_t column_index) const {
	if (column_index >= this->variant_positions.size()) {
		throw runtime_error("ColumnIndexer::get_variant_id: column index does not exist.");
	} else {
		return this->variant_positions.at(column_index);
	}
}

size_t ColumnIndexer::size() const {
	return this->variant_positions.size();
}

unsigned short ColumnIndexer::nr_paths() const {
	return this->paths.size();
}

unsigned short ColumnIndexer::get_path(unsigned short path_index) const {
	if (path_index >= this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_path: path_index does not exist.");
	} else {
		return this->paths.at(path_index);
	}
}

unsigned short ColumnIndexer::get_allele (unsigned short path_index, size_t column_index) const {
	if (path_index >= this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_allele: path_index does not exist.");
	}
	if (column_index >= this->variant_positions.size()) {
		throw runtime_error("ColumnIndex::get_allele: column_index does not exist.");
	}
	size_t idx = column_index * this->paths.size() + path_index;
	return this->allele_by_column_path.at(idx);
}

pair<unsigned short,unsigned short> ColumnIndexer::get_path_ids_at (size_t position) const {
	if (position >= this->paths.size()*this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_path_ids_at: index out of bounds.");
	}
	size_t p_id1 = position % this->nr_paths();
	size_t p_id2 = (position / this->nr_paths()) % this->nr_paths();
	return pair<unsigned short,unsigned short> (p_id2, p_id1);
}
