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
        HMMColumn* last_column = this->viterbi_columns[column_count-1];
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
                best_index = this->viterbi_backtrace_columns[column_index]->at(best_index);
                column_index -= 1;
        }
}

void HMM::compute_forward_column(size_t column_index) {
        assert(column_index < this->column_indexer->size());
        size_t variant_id = this->column_indexer->get_variant_id(column_index);
        if (this->forward_columns[column_index] != nullptr) return;

        unsigned short nr_paths = this->column_indexer->nr_paths();
        std::vector<unsigned short> col_alleles(nr_paths);
        for (unsigned short p = 0; p < nr_paths; ++p) col_alleles[p] = this->column_indexer->get_allele(p, column_index);

        HMMColumn* previous_column = nullptr;
        TransitionProbabilityComputer* transition_probability_computer = nullptr;

        if (column_index > 0) {
                previous_column = this->forward_columns[column_index-1];
                size_t prev_index = this->column_indexer->get_variant_id(column_index-1);
                size_t cur_index = this->column_indexer->get_variant_id(column_index);
                size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
                size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
                transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
        }

        HMMColumn* current_column = new HMMColumn();
        current_column->column.resize((size_t)nr_paths * nr_paths);
        current_column->row_sums.assign(nr_paths, 0.0);
        EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);

        const double* prev_buf = (column_index > 0) ? previous_column->column.data() : nullptr;
        const double* prev_row_sums = (column_index > 0) ? previous_column->row_sums.data() : nullptr;

        // Use cached row sums if available, else compute.
        std::vector<double> helper_i_vec;
        const double* helper_i = nullptr;
        if (column_index > 0) {
                if (prev_row_sums && !previous_column->row_sums.empty()) {
                        helper_i = prev_row_sums;
                } else {
                        helper_i_vec.assign(nr_paths, 0.0);
                        for (unsigned short p = 0; p < nr_paths; ++p) {
                                double row_sum = 0.0;
                                const double* row_ptr = prev_buf + (size_t)p * nr_paths;
                                for (unsigned short q = 0; q < nr_paths; ++q) row_sum += row_ptr[q];
                                helper_i_vec[p] = row_sum;
                        }
                        helper_i = helper_i_vec.data();
                }
        }

        double trans0 = 0.0, trans1 = 0.0, trans2 = 0.0;
        if (column_index > 0) {
                trans0 = transition_probability_computer->compute_transition_prob(0);
                trans1 = transition_probability_computer->compute_transition_prob(1);
                trans2 = transition_probability_computer->compute_transition_prob(2);
        }

        double A = trans0 - 2.0 * trans1 + trans2;
        double B = trans1 - trans2;
        double C = trans2; // helper_ij == 1.0

        double normalization_sum = 0.0;
        double* cur_buf = current_column->column.data();
        double* cur_row_sums = current_column->row_sums.data();

        for (unsigned short p = 0; p < nr_paths; ++p) {
                unsigned short a1 = col_alleles[p];
                double hi = (column_index > 0) ? helper_i[p] : 0.0;
                double row_sum = 0.0;
                for (unsigned short q = p; q < nr_paths; ++q) {
                        double previous_cell;
                        if (column_index > 0) {
                                double prev_pq = prev_buf[(size_t)p * nr_paths + q];
                                double hj = helper_i[q];
                                previous_cell = prev_pq * A + (hi + hj) * B + C;
                        } else {
                                previous_cell = 1.0;
                        }

                        double ep = emission_probability_computer.get_emission_probability(a1, col_alleles[q]);
                        double cell = previous_cell * ep;
                        
                        cur_buf[(size_t)p * nr_paths + q] = cell;
                        row_sum += cell;
                        if (p != q) {
                                cur_buf[(size_t)q * nr_paths + p] = cell;
                                cur_row_sums[q] += cell;
                                normalization_sum += 2.0 * cell;
                        } else {
                                normalization_sum += cell;
                        }
                }
                cur_row_sums[p] += row_sum;
        }

        if (normalization_sum > 0.0) {
                double inv_sum = 1.0 / normalization_sum;
                for (size_t j = 0; j < (size_t)nr_paths * nr_paths; ++j) cur_buf[j] *= inv_sum;
                for (unsigned short p = 0; p < nr_paths; ++p) cur_row_sums[p] *= inv_sum;
                current_column->forward_normalization_sum = normalization_sum;
        } else {
                double uniform = 1.0 / (double) (nr_paths * nr_paths);
                fill(current_column->column.begin(), current_column->column.end(), uniform);
                fill(current_column->row_sums.begin(), current_column->row_sums.end(), 1.0 / (double)nr_paths);
                current_column->forward_normalization_sum = 1.0;
        }

        this->forward_columns[column_index] = current_column;
        if (transition_probability_computer != nullptr) delete transition_probability_computer;
}

void HMM::compute_backward_column(size_t column_index) {
        size_t column_count = this->column_indexer->size();
        assert(column_index < column_count);
        size_t variant_id = this->column_indexer->get_variant_id(column_index);

        unsigned short nr_paths = column_indexer->nr_paths();
        std::vector<unsigned short> col_alleles(nr_paths);
        for (unsigned short p = 0; p < nr_paths; ++p) col_alleles[p] = this->column_indexer->get_allele(p, column_index);

        TransitionProbabilityComputer* transition_probability_computer = nullptr;
        EmissionProbabilityComputer* next_emission_computer = nullptr;
        HMMColumn* forward_column = this->forward_columns[column_index];

        if (column_index < column_count-1) {
                assert (this->previous_backward_column != nullptr);
                size_t cur_v_id = this->column_indexer->get_variant_id(column_index);
                size_t next_v_id = this->column_indexer->get_variant_id(column_index+1);
                size_t cur_pos = this->unique_kmers->at(cur_v_id)->get_variant_position();
                size_t next_pos = this->unique_kmers->at(next_v_id)->get_variant_position();
                transition_probability_computer = new TransitionProbabilityComputer(cur_pos, next_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
                next_emission_computer = new EmissionProbabilityComputer(this->unique_kmers->at(next_v_id), this->probabilities);

                if (forward_column == nullptr) {
                        size_t k = (size_t)sqrt(column_count);
                        size_t start = (column_index / k) * k;
                        for (size_t j = start + 1; j <= column_index; ++j) compute_forward_column(j);
                        forward_column = this->forward_columns[column_index];
                }
        }

        vector<double> helper_i(nr_paths, 0.0);
        const double* prev_back_buf = (column_index < column_count - 1) ? this->previous_backward_column->column.data() : nullptr;

        std::vector<unsigned short> next_alleles;
        double helper_ij = 0.0;
        if (column_index < column_count - 1) {
                next_alleles.resize(nr_paths);
                for (unsigned short p = 0; p < nr_paths; ++p) next_alleles[p] = this->column_indexer->get_allele(p, column_index + 1);

                for (unsigned short p = 0; p < nr_paths; ++p) {
                        double row_sum = 0.0;
                        unsigned short a1 = next_alleles[p];
                        for (unsigned short q = p; q < nr_paths; ++q) {
                                double term = prev_back_buf[(size_t)p * nr_paths + q] * next_emission_computer->get_emission_probability(a1, next_alleles[q]);
                                row_sum += term;
                                if (p != q) {
                                        helper_i[q] += term;
                                        helper_ij += 2.0 * term;
                                } else {
                                        helper_ij += term;
                                }
                        }
                        helper_i[p] += row_sum;
                }
        }

        HMMColumn* current_column = new HMMColumn();
        current_column->column.resize((size_t)nr_paths * nr_paths);
        current_column->row_sums.assign(nr_paths, 0.0);

        double trans0 = 0.0, trans1 = 0.0, trans2 = 0.0;
        if (column_index < column_count - 1) {
                trans0 = transition_probability_computer->compute_transition_prob(0);
                trans1 = transition_probability_computer->compute_transition_prob(1);
                trans2 = transition_probability_computer->compute_transition_prob(2);
        }

        double A = trans0 - 2.0 * trans1 + trans2;
        double B = trans1 - trans2;
        double C = helper_ij * trans2;

        double normalization_sum = 0.0;
        double normalization_f_b = 0.0;
        const double* fwd_buf = forward_column ? forward_column->column.data() : nullptr;
        double* cur_back_buf = current_column->column.data();
        double* cur_back_row_sums = current_column->row_sums.data();
        double fwd_norm = forward_column ? forward_column->forward_normalization_sum : 1.0;

        unsigned short max_allele = 0;
        for (auto a : col_alleles) if (a > max_allele) max_allele = a;
        const size_t N = (size_t)max_allele + 1;
        constexpr size_t STACK_NN = 256;
        double stack_acc[STACK_NN];
        std::vector<double> heap_acc;
        double* acc_buf;
        if (N * N <= STACK_NN) {
                std::fill(stack_acc, stack_acc + N * N, 0.0);
                acc_buf = stack_acc;
        } else {
                heap_acc.assign(N * N, 0.0);
                acc_buf = heap_acc.data();
        }

        for (unsigned short p = 0; p < nr_paths; ++p) {
                unsigned short a1 = col_alleles[p];
                double hi = helper_i[p];
                double row_sum = 0.0;
                for (unsigned short q = p; q < nr_paths; ++q) {
                        double backward_prob;
                        if (column_index < column_count - 1) {
                                double h_cell = prev_back_buf[(size_t)p * nr_paths + q] * next_emission_computer->get_emission_probability(next_alleles[p], next_alleles[q]);
                                double hj = helper_i[q];
                                backward_prob = h_cell * A + (hi + hj) * B + C;
                        } else {
                                backward_prob = 1.0;
                        }

                        cur_back_buf[(size_t)p * nr_paths + q] = backward_prob;
                        row_sum += backward_prob;
                        if (p != q) {
                                cur_back_buf[(size_t)q * nr_paths + p] = backward_prob;
                                cur_back_row_sums[q] += backward_prob;
                                normalization_sum += 2.0 * backward_prob;
                        } else {
                                normalization_sum += backward_prob;
                        }

                        if (forward_column != nullptr) {
                                double posterior = fwd_buf[(size_t)p * nr_paths + q] * backward_prob;
                                unsigned short a2 = col_alleles[q];
                                if (p != q) {
                                        acc_buf[(size_t)a1 * N + a2] += posterior * fwd_norm;
                                        acc_buf[(size_t)a2 * N + a1] += posterior * fwd_norm;
                                        normalization_f_b += 2.0 * posterior;
                                } else {
                                        acc_buf[(size_t)a1 * N + a2] += posterior * fwd_norm;
                                        normalization_f_b += posterior;
                                }
                        }
                }
                cur_back_row_sums[p] += row_sum;
        }

        {
                auto& result_for_variant = this->genotyping_result.at(variant_id);
                for (size_t a1 = 0; a1 < N; ++a1) {
                        for (size_t a2 = 0; a2 < N; ++a2) {
                                double v = acc_buf[a1 * N + a2];
                                if (v != 0.0) result_for_variant.add_to_likelihood((unsigned short)a1, (unsigned short)a2, v);
                        }
                }
        }

        if (normalization_sum > 0.0) {
                double inv_sum = 1.0 / normalization_sum;
                for (size_t j = 0; j < (size_t)nr_paths * nr_paths; ++j) cur_back_buf[j] *= inv_sum;
                for (unsigned short p = 0; p < nr_paths; ++p) cur_back_row_sums[p] *= inv_sum;
        } else {
                double uniform = 1.0 / (double) (nr_paths * nr_paths);
                fill(current_column->column.begin(), current_column->column.end(), uniform);
                fill(current_column->row_sums.begin(), current_column->row_sums.end(), 1.0 / (double)nr_paths);
        }

        if (this->previous_backward_column != nullptr) delete this->previous_backward_column;
        this->previous_backward_column = current_column;
        if (next_emission_computer != nullptr) delete next_emission_computer;

        if (this->forward_columns[column_index] != nullptr) {
                delete this->forward_columns[column_index];
                this->forward_columns[column_index] = nullptr;
        }
        if (transition_probability_computer != nullptr) delete transition_probability_computer;
}

void HMM::compute_viterbi_column(size_t column_index) {
        assert(column_index < this->column_indexer->size());
        size_t variant_id = this->column_indexer->get_variant_id(column_index);
        if (this->viterbi_columns[column_index] != nullptr) return;
        HMMColumn* previous_column = nullptr;
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
        HMMColumn* current_column = new HMMColumn();
        EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);
        double normalization_sum = 0.0;
        vector<size_t>* backtrace_column = new vector<size_t>();
        for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
                for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
                        unsigned short path1 = this->column_indexer->get_path(path_id1);
                        unsigned short path2 = this->column_indexer->get_path(path_id2);
                        double previous_cell = 0.0;
                        if (column_index > 0) {
                                size_t j = 0;
                                double max_value = 0.0;
                                size_t max_index = 0;
                                for (unsigned short prev_path_id1 = 0; prev_path_id1 < nr_paths; ++prev_path_id1) {
                                        for (unsigned short prev_path_id2 = 0; prev_path_id2 < nr_paths; ++prev_path_id2) {
                                                unsigned short prev_path1 = this->column_indexer->get_path(prev_path_id1);
                                                unsigned short prev_path2 = this->column_indexer->get_path(prev_path_id2);
                                                double prev_prob = previous_column->column.at(j);
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
                        unsigned short allele1 = this->column_indexer->get_allele(path_id1, column_index);
                        unsigned short allele2 = this->column_indexer->get_allele(path_id2, column_index);
                        double emission_prob = emission_probability_computer.get_emission_probability(allele1,allele2);
                        double current_cell = previous_cell * emission_prob;
                        current_column->column.push_back(current_cell);
                        normalization_sum += current_cell;
                }
        }
        if (normalization_sum > 0.0) {
                transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(), bind(divides<double>(), placeholders::_1, normalization_sum));
        } else {
                double uniform = 1.0 / (double) current_column->column.size();
                transform(current_column->column.begin(), current_column->column.end(), current_column->column.begin(),  [uniform](double c) -> double {return uniform;});
        }
        this->viterbi_columns[column_index] = current_column;
        this->viterbi_backtrace_columns[column_index] = backtrace_column;
        if (transition_probability_computer != nullptr) delete transition_probability_computer;
}

vector<GenotypingResult> HMM::get_genotyping_result() const {
        return this->genotyping_result;
}

vector<GenotypingResult> HMM::move_genotyping_result() {
        return move(this->genotyping_result);
}

void HMM::combine_likelihoods(HMM& other) {
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
