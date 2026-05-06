#ifndef PROBABILITYCOMPUTER_HPP
#define PROBABILITYCOMPUTER_HPP

#include <vector>
#include <stdio.h>
#include <stdlib.h>

/** 
* Computes probabilities for kmer copy numbers.
**/

class ProbabilityComputer {
public:
	ProbabilityComputer();
	ProbabilityComputer(double mean_cn0, double mean_cn1, double mean_cn2);
	void set_parameters(double mean_cn0, double mean_c1, double mean_c2);
	double get_probability (size_t cn, unsigned int value) const;
private:
	std::vector<double> means;
	double poisson(double mean, unsigned int value) const;
	double geometric(double p, unsigned int value) const;
};
#endif // PROBABILITYCOMPUTER_HPP
