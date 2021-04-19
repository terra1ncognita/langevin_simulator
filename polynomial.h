#pragma once
#include <vector>
#include <array>
#include <string>
#include "json.hpp"
#include <algorithm>

using json = nlohmann::json;

void check_key_json(const json& jsObj, std::string key);

//template<typename T, std::size_t N>
//std::ostream& operator<< (std::ostream& out, const std::array<T, N>& vec);
//
//template<typename T>
//std::ostream& operator<< (std::ostream& out, const std::vector<T>& vec);



template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& vec) {
	for (auto c : vec) {
		out << c << " ";
	}
	return out;
}


template<typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const std::array<T, N>& vec) {
	out << vec.size() << " " << vec.at(0) << std::endl;
	for (size_t i = 0; i < vec.size(); ++i) {
		out << vec[i] << " ";
	}
	//std::copy(vec.cbegin(), vec.cend(), std::ostream_iterator<T>(out, " "));
	return out;
}


class LaurentPolynomial {
public:
	LaurentPolynomial() : n{ 0 }, m{ 0 }, pole_{ 0 }, coeffs{}  {}
	LaurentPolynomial(int n, int m, double pole, std::vector<double> coeffs) : n(n), m(m), pole_(pole), coeffs(coeffs) {
		check_validity();
	}

	LaurentPolynomial(const json& jsObj, std::string key) {
		check_key_json(jsObj, key);

		n = static_cast<int>(jsObj[key]["n"]);
		m = static_cast<int>(jsObj[key]["m"]);
		pole_ = static_cast<double>(jsObj[key]["pole"]);
		coeffs = static_cast<std::vector<double>>(jsObj[key]["coeffs"]);

		check_validity();
	}

	double operator()(double x) const {
		double frac_sum = coeffs[0], res = coeffs.back(), xmp = x - pole_;

		for (int i = 1; i <= -m; ++i) {
			frac_sum = frac_sum * xmp + coeffs[i];
		}
		for (int i = n - m - 1; i >= 1 - m; --i) {
			res = res * x + coeffs[i];
		}
		return res * x + frac_sum;
	}

	std::string to_string() const {
		std::stringstream ss;
		ss << "LaurentPolynomial(n=" << n << ", m=" << m << ", pole=" << pole_ << ", coeffs=[" << coeffs << "]";
		return ss.str();
	}

	double pole() const {
		return pole_;
	}

private:
	int n, m;
	double pole_;
	std::vector<double> coeffs;

	void check_validity() {
		if (m > 0) {
			throw std::runtime_error("m should be non-positive");
		}
		if (n - m + 1 != coeffs.size()) {
			throw std::runtime_error("n, m and coeff size are not in agreement!");
		}
	}
};

std::ostream& operator<< (std::ostream& out, const LaurentPolynomial& lp);

template<typename T>
double evaluate_polynomial(const T& coeffs, double x) {
	double res = coeffs.back();
	for (int i = coeffs.size() - 2; i >= 0; --i) {
		res = res * x + coeffs[i];
	}
	return res;
}


class RationalFunction {
public:
	RationalFunction() : p{}, q{} {}
	RationalFunction(std::vector<double> pp, std::vector<double> qq) : p{ pp }, q{ qq } {}

	double operator() (double x) const {
		double num = evaluate_polynomial(p, x);
		double denom = evaluate_polynomial(q, x);
		return num / denom;
	}
private:
	std::vector<double> p, q;
};


class CubicSpline {
public:
	CubicSpline() : n{ 0 }, left_boundary{ 0.0 }, right_boundary{ 0.0 } {};

	static std::vector<double> _extract_sorted_points(const json& jsObj) {
		std::cout << "in init" << std::endl;
		std::vector<double> _points(jsObj["x"]);
		std::sort(_points.begin(), _points.end());
		return _points;
	}

	CubicSpline(const json& jsObj) :
		points{ CubicSpline::_extract_sorted_points(jsObj) },
		left_boundary{ points.front() },
		right_boundary { points.back() },
		n { points.size() } 
	{
		coef.reserve(n - 1);
		for (size_t i = 0; i < 4; ++i) {
			std::string i_str = std::to_string(i);
			if (jsObj[i_str].size() != n - 1) {
				std::string err_msg = "coeff size is not equal to n-1";
				throw std::runtime_error(err_msg);
			}
			for (int j = 0; j < n-1; ++j) {
				coef[j][i] = jsObj[i_str][j];
			}
		}
	}

	int find_domain(double x) const {
		if (x < points[0]){
			return -1;
		}
		if (x > points.back()) {
			return -2;
		}

		auto iter = std::upper_bound(points.begin(), points.end(), x);
		if (iter != points.end()) { // found
			return std::distance(points.begin(), std::prev(iter));
		}
		return -3;
	}

	double operator()(double x) const {
		int pos = find_domain(x);
		if (pos < 0) {
			std::string err_msg = "out of boundaries!";
			throw std::runtime_error(err_msg);
		}
		return evaluate_polynomial(coef[pos], x - points[pos]);
	}

	
private:
	const std::vector<double> points;
	const int n;
	std::vector<std::array<double, 4>> coef;
public:
	const double left_boundary, right_boundary;
};
