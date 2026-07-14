#include "catch.hpp"
#include "utils.hpp"
#include "../src/commands.hpp"
#include "../src/probabilitytable.hpp"
#include "../src/hmm.hpp"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cereal/archives/binary.hpp>

using namespace std;

static string replace_all(string text, const string& needle, const string& replacement) {
	size_t pos = 0;
	while ((pos = text.find(needle, pos)) != string::npos) {
		text.replace(pos, needle.size(), replacement);
		pos += replacement.size();
	}
	return text;
}

static void rewrite_fixture(const string& src, const string& dst, const string& needle, const string& replacement) {
	ifstream in(src);
	REQUIRE(in.good());
	ofstream out(dst);
	REQUIRE(out.good());

	string line;
	while (getline(in, line)) {
		out << replace_all(line, needle, replacement) << '\n';
	}
}

static vector<string> split_string(const string& text, char delim) {
	vector<string> parts;
	string part;
	istringstream iss(text);
	while (getline(iss, part, delim)) {
		parts.push_back(part);
	}
	return parts;
}

TEST_CASE("Commands run_genotype_command1", "[Commands run_genotype_command1]") {

	string precomputed_prefix = "../tests/data/index";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "../tests/data/testfull";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	double effective_N = 0.00001;
	double regularization = 0.01;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 215;
	uint64_t hash_size = 100000;
	size_t panel_size = 0;
	double recombrate = 1.26;
	bool output_panel = false;

	/** (1) produce results with command **/

	run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size, panel_size, recombrate, output_panel);

	// check if output file exists
	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}
	// parse output file
	string line;
	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf", computed_lines);

	// check if output looks as expected
	REQUIRE(computed_lines.size() == 2);


	/** (2) produce results directly from internal HMM **/

	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;

	// Reconstruct expected UniqueKmers objects from file
	ifstream is("../tests/data/region_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is( is );
	archive_is(uk);

	HMM hmm(&uk.unique_kmers["chr1"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();

	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();
	vector<vector<unsigned short>> defined = {{0,1}, {0,1,2}};
	vector<string> expected_likelihoods = {};

	for(size_t i = 0; i < 2; ++i) {
		vector<double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		ostringstream all;
		pair<int,int> genotype = genotypes[i].get_specific_likelihoods(defined[i]).get_likeliest_genotype();
		all << genotype.first << "/" << genotype.second << ":";
		all << genotypes[i].get_specific_likelihoods(defined[i]).get_genotype_quality(genotype.first, genotype.second) << ":";
		all << setprecision(4) << log10(likelihoods[0]);
		for (size_t j = 1; j < likelihoods.size(); ++j) {
			all << "," << setprecision(4) << log10(likelihoods[j]);
		}
		all << ":" << uk.unique_kmers["chr1"][i]->get_coverage();
		expected_likelihoods.push_back(all.str());
	}


	/** (3) Check if results are identical **/

	for (size_t i = 0; i < expected_likelihoods.size(); ++i) {
		REQUIRE(expected_likelihoods[i] == computed_lines[i][9]);
	}
}

TEST_CASE("Commands run_single_command_chrX", "[Commands run_single_command_chrX]") {
	string chrX_vcf = "chrX-region.vcf";
	string chrX_ref = "chrX-region.fa";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "chrX-single";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	double effective_N = 0.00001;
	double regularization = 0.01;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	bool add_reference = false;
	size_t sampling_size = 215;
	uint64_t hash_size = 100000;
	size_t panel_size = 0;
	double recombrate = 1.26;
	bool output_panel = false;
	double sampling_effective_N = 0.01;
	unsigned short allele_penalty = 5;
	bool serialize_output = false;

	rewrite_fixture("../tests/data/region.vcf", chrX_vcf, "chr1", "chrX");
	rewrite_fixture("../tests/data/region.fa", chrX_ref, "chr1", "chrX");

	run_single_command(outname, readfile, chrX_ref, chrX_vcf, 31, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, add_reference, sampling_size, hash_size, panel_size, recombrate, output_panel, sampling_effective_N, allele_penalty, serialize_output);

	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}

	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf", computed_lines);
	REQUIRE(computed_lines.size() == 2);
	REQUIRE(computed_lines[0][0] == "chrX");
	REQUIRE(computed_lines[1][0] == "chrX");

	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;
	ifstream is("../tests/data/region_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is(is);
	archive_is(uk);
	uk.unique_kmers["chrX"] = uk.unique_kmers["chr1"];

	HMM hmm(&uk.unique_kmers["chrX"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();
	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();

	vector<vector<unsigned short>> defined = {{0,1}, {0,1,2}};
	vector<vector<double>> expected_gls;
	vector<string> expected_gt;
	vector<string> expected_gq;
	vector<string> expected_kc;
	for (size_t i = 0; i < 2; ++i) {
		if (genotypes[i].contains_no_likelihoods()) {
			genotypes[i].add_to_likelihood(0, 0, 1.0);
		}
		vector<double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		pair<int,int> genotype = genotypes[i].get_specific_likelihoods(defined[i]).get_likeliest_genotype();
		if ((genotype.first != -1) && (genotype.second != -1)) {
			expected_gt.push_back(to_string(genotype.first) + "/" + to_string(genotype.second));
			expected_gq.push_back(to_string(genotypes[i].get_specific_likelihoods(defined[i]).get_genotype_quality(genotype.first, genotype.second)));
		} else {
			expected_gt.push_back(".");
			expected_gq.push_back(".");
		}
		vector<double> gl_values;
		gl_values.push_back(log10(likelihoods[0]));
		for (size_t j = 1; j < likelihoods.size(); ++j) {
			gl_values.push_back(log10(likelihoods[j]));
		}
		expected_gls.push_back(gl_values);
		expected_kc.push_back(to_string(uk.unique_kmers["chrX"][i]->get_coverage()));
	}

	for (size_t i = 0; i < expected_gls.size(); ++i) {
		vector<string> fields = split_string(computed_lines[i][9], ':');
		REQUIRE(fields.size() == 4);
		REQUIRE(fields[0] == expected_gt[i]);
		REQUIRE(fields[1] == expected_gq[i]);
		REQUIRE(fields[3] == expected_kc[i]);

		vector<string> gl_tokens = split_string(fields[2], ',');
		REQUIRE(gl_tokens.size() == expected_gls[i].size());
		for (size_t j = 0; j < gl_tokens.size(); ++j) {
			REQUIRE(fabs(stod(gl_tokens[j]) - expected_gls[i][j]) < 0.1);
		}
	}
}

TEST_CASE("Commands run_genotype_command2", "[Commands run_genotype_command2]") {
	string precomputed_prefix = "../tests/data/index";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "../tests/data/testsampled";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	double effective_N = 0.00001;
	double regularization = 0.01;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 0;
	uint64_t hash_size = 100000;
	size_t panel_size = 5;
	double recombrate = 1.26;
	bool output_panel = false;

	/** (1) produce results with command **/

	run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size, panel_size, recombrate, output_panel);
	// check if output file exists
	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}
	// parse output file
	string line;
	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf", computed_lines);

	// check if output looks as expected
	REQUIRE(computed_lines.size() == 2);


	/** (2) produce results directly from internal HMM **/

	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;

	// Reconstruct expected UniqueKmers objects from file
	ifstream is("../tests/data/region2_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is( is );
	archive_is(uk);

	HMM hmm(&uk.unique_kmers["chr1"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();

	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();
	vector<vector<unsigned short>> defined = {{0,1}, {0,1,2}};
	vector<string> expected_likelihoods = {};

	for(size_t i = 0; i < 2; ++i) {
		if (genotypes[i].contains_no_likelihoods()) {
			genotypes[i].add_to_likelihood(0,0,1.0);
		}
		vector<double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		ostringstream all;
		pair<int,int> genotype = genotypes[i].get_specific_likelihoods(defined[i]).get_likeliest_genotype();
		if ((genotype.first != -1) && (genotype.second != -1)) {
			all << genotype.first << "/" << genotype.second << ":";
			all << genotypes[i].get_specific_likelihoods(defined[i]).get_genotype_quality(genotype.first, genotype.second) << ":";
		} else {
			all << ".:.:";
		}
		all << setprecision(4) << log10(likelihoods[0]);
		for (size_t j = 1; j < likelihoods.size(); ++j) {
			all << "," << setprecision(4) << log10(likelihoods[j]);
		}

		all << ":" << uk.unique_kmers["chr1"][i]->get_coverage();
		expected_likelihoods.push_back(all.str());
	}


	/** (3) Check if results are identical **/

	for (size_t i = 0; i < expected_likelihoods.size(); ++i) {
		REQUIRE(expected_likelihoods[i] == computed_lines[i][9]);
	}
}
