#include <iostream>
#include <fstream>
#include "fastareader.hpp"

using namespace std;

FastaReader::FastaReader(const string& filename) {
	parse_file(filename);
	cerr << "Found " << this->name_to_sequence.size() << " chromosome(s) from the reference file." << endl;
}


// std::map<std::string, DnaSequence*> name_to_sequence;
void FastaReader::parse_file(const string& filename) {
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("FastaReader::parse_file: reference file cannot be opened.");
	}
	string line;
	shared_ptr<DnaSequence> dna_seq = nullptr;
	while (getline(file, line)) {
		if (line.size() == 0) continue;
		size_t start = line.find_first_not_of(" \t\r\n");
		if (start == string::npos) continue;
		size_t end = line.find_last_not_of(" \t\r\n");
		if (end + 1 < line.size()) line.erase(end + 1);
		if (start > 0) line.erase(0, start);
		if (line.size() == 0) continue;
		// sequence name ?
		if (line[0] == '>') {
			start = line.find_first_not_of(" \t", 1);
			end = line.find_first_of(" \t", start);
			if (end == string::npos) end = line.size();
			string name = line.substr(start,end-start);

			dna_seq = make_shared<DnaSequence>();
			this->name_to_sequence[name] = dna_seq;
		} else {
			if (dna_seq == nullptr) {
				throw runtime_error("FastaReader::parse_file: file is malformatted.");
			} else {
				dna_seq->append(line);
			}
		}
	}
}


bool FastaReader::contains_name(const string& name) const {
	auto it = this->name_to_sequence.find(name);
	return (it != this->name_to_sequence.end());
}

size_t FastaReader::get_size_of(const string& name) const {
	auto it = this->name_to_sequence.find(name);
	if (it != this->name_to_sequence.end()) return it->second->size();
	throw runtime_error("FastaReader::get_size_of: chromosome " + name + " is not present in FASTA-file.");
}

void FastaReader::get_sequence_names(vector<string>& names) const {
	for (auto it = this->name_to_sequence.begin(); it != this->name_to_sequence.end(); ++it) {
		names.push_back(it->first);
	}
}

size_t FastaReader::get_total_kmers(size_t kmer_size) const {
	size_t total_kmers = 0;
	for (auto it = this->name_to_sequence.begin(); it != this->name_to_sequence.end(); ++it) {
		total_kmers += it->second->size() - kmer_size + 1;
	}
	return total_kmers;
}

void FastaReader::get_subsequence(const string& name, size_t start, size_t end, string& result) const {
	auto it = this->name_to_sequence.find(name);
	if (it != this->name_to_sequence.end()) {
		it->second->substr(start, end, result);
		return;
	}
	throw runtime_error("FastaReader::get_subsequence (string): chromosome " + name + " is not present in FASTA-file.");
}

void FastaReader::get_subsequence(const string& name, size_t start, size_t end, DnaSequence& result) const {
	auto it = this->name_to_sequence.find(name);
	if (it != this->name_to_sequence.end()) {
		it->second->substr(start, end, result);
		return;
	}
	throw runtime_error("FastaReader::get_subsequence (DnaSequence): chromosome " + name + " is not present in FASTA-file.");
}

FastaReader FastaReader::extract_name(const string& name) {
	auto it = this->name_to_sequence.find(name);
	if (it != this->name_to_sequence.end()) {
		// sequence name exists. Extract it
		FastaReader extracted;
		extracted.name_to_sequence[name] = move(it->second);
		this->name_to_sequence.erase(it);
		return extracted;
	} else {
		throw runtime_error("FastaReader::extract_name: chromosome " + name + " is not present in FASTA-file.");		
	}
}
