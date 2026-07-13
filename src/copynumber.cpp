#include <stdexcept>
#include <sstream>
#include "copynumber.hpp"

using namespace std;

CopyNumber::CopyNumber()
:probabilities({1.0, 0.0, 0.0})
{}

CopyNumber::CopyNumber(double cn_0, double cn_1, double cn_2)
	:probabilities({cn_0, cn_1, cn_2})
{}

CopyNumber::CopyNumber(double cn_0, double cn_1, double cn_2, double regularization_const)
	:probabilities({0.0, 0.0, 0.0})
{
	double sum = cn_0 + cn_1 + cn_2 + 3.0 * regularization_const;
	this->probabilities[0] = (cn_0 + regularization_const) / sum;
	this->probabilities[1] = (cn_1 + regularization_const) / sum;
	// Keep the original implementation's derivation and operation order.
	this->probabilities[2] = 1.0 - probabilities[0] - probabilities[1];
}

double CopyNumber::get_probability_of(int cn) const {
	if( (cn < 0) || (cn > 2) ){
		ostringstream oss;
		oss << "CopyNumber::get_probability_of: Invalid copy number: " << cn;
		throw runtime_error(oss.str());
	}
	return probabilities[cn];
}

bool CopyNumber::operator==(const CopyNumber &other) const{
	return this->probabilities == other.probabilities;
}

bool CopyNumber::operator!=(const CopyNumber &other) const {
	return !(*this == other);
}
