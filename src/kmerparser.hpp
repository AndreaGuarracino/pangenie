#include <vector>
#include <string>
#include <string_view>
#include <memory>
#include "kmercounter.hpp"

void parse(std::vector<std::string>& result, const std::string& line, char sep);

void parse_kmer_line(const std::string& line, std::string& chromosome, size_t& start, std::vector<std::string>& kmers, std::vector<std::string>& flanking_kmers, bool& is_header);

void parse_kmer_line(std::string_view line, std::string_view& chromosome, size_t& start, std::vector<std::string_view>& kmers, std::vector<std::string_view>& flanking_kmers, bool& is_header);

unsigned short compute_local_coverage(std::span<const size_t> counts, size_t kmer_coverage);
