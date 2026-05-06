#include <utility>
#include <math.h>
#include <cassert>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "hmm.hpp"
#include "emissionprobabilitycomputer.hpp"

#include <iostream>

using namespace std;


void print_column(vector<double>* column, ColumnIndexer* indexer) {
	for (size_t i = 0; i < column->size(); ++i) {
		pair<size_t,size_t> paths = indexer->get_path_ids_at(i);
		cout << setprecision(15) << column->at(i) << " paths: " << paths.first << " " <<  paths.second << endl;
	}
	cout << "" << endl;
}


HMM::HMM(vector<shared_ptr<UniqueKmers>>* unique_kmers, ProbabilityTable* probabilities, bool run_genotyping, bool run_phasing, double recombrate, bool uniform, double effective_N, vector<unsigned short>* only_paths, bool normalize)
	: unique_kmers(unique_kmers),
	 probabilities(probabilities),
	 genotyping_result(unique_kmers->size()),
	 recombrate(recombrate),
	 uniform(uniform),
	 effective_N(effective_N)
{
	this->column_indexer = new ColumnIndexer(unique_kmers, only_paths);
	this->previous_backward_column = nullptr;

	if (run_genotyping) {
		compute_forward_prob();
		compute_backward_prob();

		if (normalize) {
			for (size_t i = 0; i < this->genotyping_result.size(); ++i) {
				genotyping_result[i].normalize();
			}
		}
	}

	if (run_phasing) {
		compute_viterbi_path();
	}

	// delete objects no longer needed to save space
	init(this->forward_columns,0);
	if (this->previous_backward_column != nullptr){
		delete this->previous_backward_column;
		this->previous_backward_column = nullptr;
	}
	init(this->viterbi_columns,0);
	init(this->viterbi_backtrace_columns,0);
	if (this->column_indexer != nullptr) {
		delete this->column_indexer;
		this->column_indexer = nullptr;
	}
}

HMM::~HMM(){
	init(this->forward_columns,0);
	if (this->previous_backward_column != nullptr) delete this->previous_backward_column;
	init(this->viterbi_columns,0);
	init(this->viterbi_backtrace_columns,0);
	if (this->column_indexer != nullptr) {
		delete this->column_indexer;
		this->column_indexer = nullptr;
	}
}

void HMM::compute_forward_prob() {
	size_t column_count = this->column_indexer->size();
	init(this->forward_columns, column_count);
	
	// forward pass
	size_t k = (size_t) sqrt(column_count);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {;
		compute_forward_column(column_index);
		// sparse table: check whether to delete previous column
		if ( (k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->forward_columns[column_index-1];
			this->forward_columns[column_index-1] = nullptr;
		}
	}
}

void HMM::compute_backward_prob() {
	size_t column_count = this->column_indexer->size();
	if (column_count == 0) return;
	if (this->previous_backward_column != nullptr) {
		delete this->previous_backward_column;
		this->previous_backward_column = nullptr;
	}

	// backward pass
	for (int column_index = column_count-1; column_index >= 0; --column_index) {
		compute_backward_column(column_index);
	}
	
	// store the number of unique kmers and coverage
	for (size_t i = 0; i < this->unique_kmers->size(); ++i) {
		this->genotyping_result.at(i).set_unique_kmers(this->unique_kmers->at(i)->size());
		this->genotyping_result.at(i).set_coverage(this->unique_kmers->at(i)->get_coverage());
	}
}

void HMM::compute_viterbi_path() {
	size_t column_count = this->column_indexer->size();
	if (column_count == 0) return;
	init(this->viterbi_columns, column_count);
	init(this->viterbi_backtrace_columns, column_count);

	// perform viterbi algorithm
	size_t k = (size_t) sqrt(column_count);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {
		compute_viterbi_column(column_index);
		// sparse table: check whether to delete previous column
		if ((k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->viterbi_columns[column_index-1];
			this->viterbi_columns[column_index-1] = nullptr;
			delete this->viterbi_backtrace_columns[column_index-1];
			this->viterbi_backtrace_columns[column_index-1] = nullptr;
		}
	}

	// find best value (+ index) in last column
	size_t best_index = 0;
	double best_value = 0.0;
	HMMColumn* last_column = this->viterbi_columns.at(column_count-1);
	assert (last_column != nullptr);
	for (size_t i = 0; i < last_column->column.size(); ++i) {
		double entry = last_column->column.at(i);
		if (entry >= best_value) {
			best_value = entry;
			best_index = i;
		}
	}

	// backtracking
	size_t column_index = column_count - 1;
	while (true) {
		pair<unsigned short, unsigned short> path_ids = this->column_indexer->get_path_ids_at(best_index);
		unsigned short allele1 = this->column_indexer->get_allele (path_ids.first, column_index);
		unsigned short allele2 = this->column_indexer->get_allele (path_ids.second, column_index);

		// columns might have to be re-computed
		if (this->viterbi_backtrace_columns[column_index] == nullptr) {
			size_t j = column_index / k*k;
			assert (this->viterbi_columns[j] != nullptr);
			for (j = j+1; j<=column_index; ++j) {
				compute_viterbi_column(j);
			}
		}

		// store resulting haplotypes
		size_t variant_id = this->column_indexer->get_variant_id(column_index);
		this->genotyping_result.at(variant_id).add_first_haplotype_allele(allele1);
		this->genotyping_result.at(variant_id).add_second_haplotype_allele(allele2);
		this->genotyping_result.at(column_index).set_unique_kmers(this->unique_kmers->at(column_index)->size());
		this->genotyping_result.at(column_index).set_coverage(this->unique_kmers->at(column_index)->get_coverage());

		if (column_index == 0) break;

		// update best index 
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);
		column_index -= 1;
	}
}

void HMM::compute_forward_column(size_t column_index) {
	// NOTE: this implementation assumes that all variant positions are covered by the same set of paths
	assert(column_index < this->column_indexer->size());
	size_t variant_id = this->column_indexer->get_variant_id(column_index);

	// check whether column was computed already
	if (this->forward_columns[column_index] != nullptr) return;

	// get previous column and previous path ids (if existent)
	HMMColumn* previous_column = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	
	// nr of paths
	unsigned short nr_paths = this->column_indexer->nr_paths();
	
	if (column_index > 0) {
		previous_column = this->forward_columns[column_index-1];
		size_t prev_index = this->column_indexer->get_variant_id(column_index-1);
		size_t cur_index = this->column_indexer->get_variant_id(column_index);
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
	}

	// construct new column
	HMMColumn* current_column = new HMMColumn();
	// stage-1 microopt: pre-size to avoid push_back capacity checks in hot loop.
	current_column->column.resize((size_t)nr_paths * nr_paths);

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);

	vector<double> helper_i(nr_paths);
	vector<double> helper_j(nr_paths);
	double helper_ij = 0.0;

	// stage-1 microopt: cache the underlying buffer to drop one indirection per iteration.
	const double* prev_buf = (column_index > 0) ? previous_column->column.data() : nullptr;

	if (column_index > 0) {
		size_t i = 0;
		for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
			for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
				double prev_forward = prev_buf[i];
				helper_i[path_id1] += prev_forward;
				helper_j[path_id2] += prev_forward;
				helper_ij += prev_forward;
				i += 1;
			}
		}
	}

	// stage-1 microopt: precompute the 3 transition probabilities once per column
	// instead of calling compute_transition_prob() 3*nr_paths^2 times.
	double trans0 = 0.0, trans1 = 0.0, trans2 = 0.0;
	if (column_index > 0) {
		trans0 = transition_probability_computer->compute_transition_prob(0);
		trans1 = transition_probability_computer->compute_transition_prob(1);
		trans2 = transition_probability_computer->compute_transition_prob(2);
	}

	// normalization
	double normalization_sum = 0.0;

	// stage-1 microopt: write directly into pre-sized buffer; no push_back in hot loop.
	double* cur_buf = current_column->column.data();

	// state index
	size_t i = 0;
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		// stage-1 microopt: hoist allele1 lookup; only depends on path_id1.
		unsigned short allele1 = this->column_indexer->get_allele(path_id1, column_index);
		double helper_i_p1 = helper_i[path_id1];
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			double previous_cell = 0.0;
			if (column_index > 0) {
				// stage-1 microopt: load prev_buf[i] once (was 3 indirections + 3 .at() bounds checks).
				double prev_i = prev_buf[i];
				double helper_j_p2 = helper_j[path_id2];
				previous_cell = trans0 * prev_i +
								trans1 * (helper_i_p1 + helper_j_p2 - 2*prev_i) +
								trans2 * (helper_ij - helper_i_p1 - helper_j_p2 + prev_i);
			} else {
				previous_cell = 1.0;
			}

			// determine alleles current paths (ids) correspond to
			unsigned short allele2 = this->column_indexer->get_allele(path_id2, column_index);
			// determine emission probability
			double emission_prob = emission_probability_computer.get_emission_probability(allele1,allele2);

			// set entry of current column
			double current_cell = previous_cell * emission_prob;
			cur_buf[i] = current_cell;
			normalization_sum += current_cell;
			i += 1;
		}
	}

	if (normalization_sum > 0.0) {
		// normalize the entries in current column to sum up to 1
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(), bind(divides<double>(), placeholders::_1, normalization_sum));
	} else {
		double uniform = 1.0 / (double) current_column->column.size();
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(),  [uniform](double c) -> double {return uniform;});
//		cerr << "Underflow in Forward pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

	// store the column
	if (normalization_sum > 0.0) {
		current_column->forward_normalization_sum = normalization_sum;
	} else {
		current_column->forward_normalization_sum = 1.0;
	}
	this->forward_columns.at(column_index) = current_column;

	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}

void HMM::compute_backward_column(size_t column_index) {
	size_t column_count = this->column_indexer->size();
	assert(column_index < column_count);
	size_t variant_id = this->column_indexer->get_variant_id(column_index);

	// get previous probabilitycomputers
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	EmissionProbabilityComputer* emission_probability_computer = nullptr;
	HMMColumn* forward_column = this->forward_columns.at(column_index);
	
	// nr of paths
	unsigned short nr_paths = column_indexer->nr_paths();

	if (column_index < column_count-1) {
		assert (this->previous_backward_column != nullptr);
		size_t prev_index = this->column_indexer->get_variant_id(column_index);
		size_t cur_index = this->column_indexer->get_variant_id(column_index+1);
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);	
		emission_probability_computer = new EmissionProbabilityComputer(this->unique_kmers->at(this->column_indexer->get_variant_id(column_index+1)), this->probabilities);

		// get forward probabilities (needed for computing posteriors
		if (forward_column == nullptr) {
			// compute index of last column stored
			size_t k = (size_t)sqrt(column_count);
			size_t next = min((size_t) ( (column_index / k) * k ), column_count-1);
			for (size_t j = next+1; j <= column_index; ++j) {
				compute_forward_column(j);
			}
		}

		forward_column = this->forward_columns.at(column_index);
		assert (forward_column != nullptr);
	}

	vector<double> helper_i(nr_paths);
	vector<double> helper_j(nr_paths);
	double helper_ij = 0.0;

	// stage-1 microopt: cache buffer pointer to drop indirection.
	const double* prev_back_buf = (column_index < column_count - 1)
		? this->previous_backward_column->column.data() : nullptr;

	if (column_index < column_count - 1) {
		size_t i = 0;
		for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
			unsigned short prev_allele1 = this->column_indexer->get_allele(path_id1, column_index + 1);
			for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
				unsigned short prev_allele2 = this->column_indexer->get_allele(path_id2, column_index + 1);
				// stage-1 microopt: compute the (prev_back * emission) product once
				// instead of three times per cell.
				double term = prev_back_buf[i] *
					emission_probability_computer->get_emission_probability(prev_allele1, prev_allele2);
				helper_i[path_id1] += term;
				helper_j[path_id2] += term;
				helper_ij += term;
				i += 1;
			}
		}
	}

	// construct new column
	HMMColumn* current_column = new HMMColumn;
	// stage-1 microopt: pre-size to avoid push_back capacity checks.
	current_column->column.resize((size_t)nr_paths * nr_paths);

	// stage-1 microopt: precompute 3 transition probabilities once.
	double trans0 = 0.0, trans1 = 0.0, trans2 = 0.0;
	if (column_index < column_count - 1) {
		trans0 = transition_probability_computer->compute_transition_prob(0);
		trans1 = transition_probability_computer->compute_transition_prob(1);
		trans2 = transition_probability_computer->compute_transition_prob(2);
	}

	// normalization
	double normalization_sum = 0.0;

	// normalization of forward-backward
	double normalization_f_b = 0.0;

	// stage-1 microopt: cache forward_column buffer pointer too.
	const double* fwd_buf = forward_column ? forward_column->column.data() : nullptr;
	double* cur_back_buf = current_column->column.data();
	double fwd_norm = forward_column ? forward_column->forward_normalization_sum : 1.0;

	// state index
	size_t i = 0;
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		// stage-1 microopt: hoist allele1 + prev_allele1 (only depend on path_id1).
		unsigned short allele1 = this->column_indexer->get_allele(path_id1, column_index);
		unsigned short prev_allele1 = (column_index < column_count - 1)
			? this->column_indexer->get_allele(path_id1, column_index + 1) : 0;
		double helper_i_p1 = helper_i[path_id1];
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get alleles on current paths
			unsigned short allele2 = this->column_indexer->get_allele(path_id2, column_index);
			double current_cell = 0.0;
			if (column_index < column_count - 1) {
				// get alleles on previous paths, assuming indexes are same as current column
				unsigned short prev_allele2 = this->column_indexer->get_allele(path_id2, column_index + 1);
				double helper_cell = prev_back_buf[i] *
					emission_probability_computer->get_emission_probability(prev_allele1, prev_allele2);
				double helper_j_p2 = helper_j[path_id2];
				// iterate over previous column (ahead of this)
				current_cell =	trans0 * helper_cell +
								trans1 * (helper_i_p1 + helper_j_p2 - 2*helper_cell) +
								trans2 * (helper_ij - helper_i_p1 - helper_j_p2 + helper_cell);
			} else {
				current_cell = 1.0;
			}
			// store computed backward prob in column
			cur_back_buf[i] = current_cell;
			normalization_sum += current_cell;

			// compute forward_prob * backward_prob
			double forward_backward_prob = fwd_buf[i] * current_cell;
			normalization_f_b += forward_backward_prob;

			// update genotype likelihood
			this->genotyping_result.at(variant_id).add_to_likelihood(allele1, allele2, forward_backward_prob * fwd_norm);
			i += 1;
		}
	}


	if (normalization_sum > 0.0) {
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(), bind(divides<double>(), placeholders::_1, normalization_sum));
	} else {
		double uniform = 1.0 / (double) current_column->column.size();
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(), [uniform](double c) -> double {return uniform;});
//		cerr << "Underflow in Backward pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

//	cout << "FORWARD COLUMN: " << endl;
//	print_column(forward_column, column_indexer);

//	cout << "BACKWARD COLUMN: "  << endl;
//	print_column(current_column, column_indexer);

	// store computed column (needed for next step)
	if (this->previous_backward_column != nullptr) {
		delete this->previous_backward_column;
		this->previous_backward_column = nullptr;
	}
	this->previous_backward_column = current_column;
	if (emission_probability_computer != nullptr) delete emission_probability_computer;

	// delete forward column as it's not needed any more
	if (this->forward_columns.at(column_index) != nullptr) {
		delete this->forward_columns.at(column_index);
		this->forward_columns.at(column_index) = nullptr;
	}

	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}


void HMM::compute_viterbi_column(size_t column_index) {
	assert(column_index < this->column_indexer->size());
	size_t variant_id = this->column_indexer->get_variant_id(column_index);

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column
	HMMColumn* previous_column = nullptr;

	// nr of paths
	unsigned short nr_paths = column_indexer->nr_paths();
	
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index-1];
		size_t prev_index = this->column_indexer->get_variant_id(column_index-1);
		size_t cur_index = this->column_indexer->get_variant_id(column_index);
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
	}

	// construct new column
	HMMColumn* current_column = new HMMColumn();

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);

	// normalization 
	double normalization_sum = 0.0;

	// backtrace table
	vector<size_t>* backtrace_column = new vector<size_t>();

	// state index
	size_t i = 0;
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			unsigned short path1 = this->column_indexer->get_path(path_id1);
			unsigned short path2 = this->column_indexer->get_path(path_id2);
			double previous_cell = 0.0;
			if (column_index > 0) {
				// previous state index
				size_t j = 0;
				double max_value = 0.0;
				size_t max_index = 0;
				// iterate over all pairs of previous paths
				for (unsigned short prev_path_id1 = 0; prev_path_id1 < nr_paths; ++prev_path_id1) {
					for (unsigned short prev_path_id2 = 0; prev_path_id2 < nr_paths; ++prev_path_id2) {
						// paths corresponding to path indices
						unsigned short prev_path1 = this->column_indexer->get_path(prev_path_id1);
						unsigned short prev_path2 = this->column_indexer->get_path(prev_path_id2);
						// probability of previous cell
						double prev_prob = previous_column->column.at(j);
						// determine transition probability
						double transition_prob = transition_probability_computer->compute_transition_prob(prev_path1, prev_path2, path1, path2);
						prev_prob *= transition_prob;
						if (prev_prob >= max_value) {
							max_value = prev_prob;
							max_index = j;
						}
						j += 1;
					}
				}
				previous_cell = max_value;
				backtrace_column->push_back(max_index);
			} else {
				previous_cell = 1.0;
			}

			// determine alleles current paths (ids) correspond to
			unsigned short allele1 = this->column_indexer->get_allele(path_id1, column_index);
			unsigned short allele2 = this->column_indexer->get_allele(path_id2, column_index);
			// determine emission probability
			double emission_prob = emission_probability_computer.get_emission_probability(allele1,allele2);
			// set entry of current column
			double current_cell = previous_cell * emission_prob;
			current_column->column.push_back(current_cell);
			normalization_sum += current_cell;
			i += 1;
		}
	}

	if (normalization_sum > 0.0) {
		// normalize the entries in current column to sum up to 1 
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(), bind(divides<double>(), placeholders::_1, normalization_sum));
	} else {
		double uniform = 1.0 / (double) current_column->column.size();
		transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(),  [uniform](double c) -> double {return uniform;});
//		cerr << "Underflow in Viterbi pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

	// store the column
	this->viterbi_columns.at(column_index) = current_column;
	if (column_index > 0) assert(backtrace_column->size() == this->column_indexer->nr_paths()*this->column_indexer->nr_paths());
	this->viterbi_backtrace_columns.at(column_index) = backtrace_column;
	
	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}

vector<GenotypingResult> HMM::get_genotyping_result() const {
	return this->genotyping_result;
}

vector<GenotypingResult> HMM::move_genotyping_result() {
	return move(this->genotyping_result);
}

void HMM::combine_likelihoods(HMM& other) {
	// TODO: implement this.
	if (this->genotyping_result.size() != other.genotyping_result.size()) {
		throw runtime_error("HMM::combine_likelihoods: HMMs to be combined must be of the same size.");
	}
	size_t index = 0;
	for (auto likelihoods : other.genotyping_result) {
		this->genotyping_result.at(index).combine(likelihoods);
		index += 1;
	}
}

void HMM::normalize() {
	for (size_t i = 0; i < this->genotyping_result.size(); ++i) {
		this->genotyping_result[i].normalize();
	}
}
