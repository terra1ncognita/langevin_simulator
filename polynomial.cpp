#include <iostream>
#include "polynomial.h"

std::ostream& operator<< (std::ostream& out, const LaurentPolynomial& lp) {
	out << lp.to_string();
	return out;
}


void check_key_json(const json& jsObj, std::string key) {
	if (jsObj[key].empty()) {
		throw std::runtime_error("Provided key is not in json");
	}
}