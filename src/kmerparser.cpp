#include<iostream>
#include<cassert>
#include <charconv>
#include <sstream>
#include "kmerparser.hpp"

using namespace std;

void parse(vector<string>& result, const string& line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

void parse_kmer_line(const string& line, string& chrom, size_t& start, vector<string>& kmers, vector<string>& flanking_kmers, bool& is_header) {
	vector<string> tokens;
	parse(tokens, line, '\t');
	assert (tokens.size() == 5);
	if (tokens[0][0] == '#') {
		is_header = true;
		return;
	}
	chrom = tokens[0];
	start = atoi(tokens[1].c_str());
	if (tokens[3] != "nan") parse(kmers, tokens[3], ',');
	if (tokens[4] != "nan") parse(flanking_kmers, tokens[4], ',');
}

namespace {

void split_views(string_view input, char separator, vector<string_view>& result) {
	size_t begin = 0;
	while (begin <= input.size()) {
		size_t end = input.find(separator, begin);
		if (end == string_view::npos) end = input.size();
		result.emplace_back(input.data() + begin, end - begin);
		if (end == input.size()) break;
		begin = end + 1;
	}
}

} // namespace

void parse_kmer_line(string_view line, string_view& chrom, size_t& start, vector<string_view>& kmers, vector<string_view>& flanking_kmers, bool& is_header) {
	string_view fields[5];
	size_t begin = 0;
	for (size_t field = 0; field < 5; ++field) {
		size_t end = line.find('\t', begin);
		if (end == string_view::npos) end = line.size();
		fields[field] = line.substr(begin, end - begin);
		begin = end + 1;
	}
	assert(begin == line.size() + 1);

	is_header = !fields[0].empty() && fields[0].front() == '#';
	if (is_header) return;

	chrom = fields[0];
	auto parsed = from_chars(fields[1].data(), fields[1].data() + fields[1].size(), start);
	if (parsed.ec != errc() || parsed.ptr != fields[1].data() + fields[1].size()) {
		throw runtime_error("parse_kmer_line: invalid variant start position.");
	}
	if (fields[3] != "nan") split_views(fields[3], ',', kmers);
	if (fields[4] != "nan") split_views(fields[4], ',', flanking_kmers);
}

unsigned short compute_local_coverage(span<const size_t> counts, size_t kmer_coverage) {
	size_t total_coverage = 0;
	size_t total_kmers = 0;
	size_t min_cov = kmer_coverage / 4;
	size_t max_cov = kmer_coverage * 4;

	for (size_t read_count : counts) {
		// ignore too extreme counts
		if ( (read_count < min_cov) || (read_count > max_cov) ) continue;
		total_coverage += read_count;
		total_kmers += 1;		
	}

	if ((total_kmers > 0) && (total_coverage > 0)){
		return total_coverage / total_kmers;
	} else {
		return kmer_coverage;
	}		
}
