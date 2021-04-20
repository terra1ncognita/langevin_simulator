#pragma once
#include <vector>
#include <array>
#include <string>
#include "json.hpp"
#include <algorithm>

using json = nlohmann::json;

void check_key_json(const json& jsObj, std::string key);

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& vec) {
	for (auto c : vec) {
		out << c << " ";
	}
	return out;
}


template<typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const std::array<T, N>& vec) {
	for (size_t i = 0; i < vec.size(); ++i) {
		out << vec[i] << " ";
	}
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
		coeffs = jsObj[key]["coeffs"].get<std::vector<double>>();

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
		std::vector<double> _points = jsObj["x"].get<std::vector<double>>();
		std::sort(_points.begin(), _points.end());
		return _points;
	}

	CubicSpline(const std::vector<double> points, const int n, std::vector<std::array<double, 4>> _coef, const double left_boundary, const double right_boundary
	) : points{ points }, n{ n }, left_boundary{ left_boundary }, right_boundary{ right_boundary }, coef{_coef} {}


	CubicSpline(const json& jsObj) :
		points{ CubicSpline::_extract_sorted_points(jsObj) },
		left_boundary{ points.front() },
		right_boundary { points.back() },
		n { points.size() } 
	{
		coef.resize(n - 1);
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

	std::size_t find_domain(double x) const {
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
		std::size_t pos = find_domain(x);
		if (pos < 0) {
			std::string err_msg = "out of boundaries!";
			throw std::runtime_error(err_msg);
		}
		return evaluate_polynomial(coef[pos], x - points[pos]);
	}

	std::ostream& print_points(std::ostream &out) const {
		out << points;
		return out;
	}

	std::ostream& print_coef(std::ostream &out) const {
		for (const auto & arr : coef) {
			out << arr << std::endl;
		}
		return out;
	}


	//
	//CubicSpline & operator=(const CubicSpline& other){
	//	std::cout << "in copy assignment" << std::endl;

	//	if (this == &other) {
	//		return *this;
	//	}

	//	std::cout << "Other coef size " << other.coef.size() << std::endl;
	//	std::cout << "Other coef " << other.coef << std::endl;

	//	this->points = other.points;
	//	this->n = other.n;
	//	this->left_boundary = other.left_boundary;
	//	this->right_boundary = other.right_boundary;

	//	coef.reserve(other.coef.size());

	//	for (auto i = 0; i < other.coef.size(); i++) {
	//		//std::unique_ptr<std::array<double, 4>> task = std::make_unique<std::array<double, 4>>(arr);
	//		//coef.push_back(*std::move(task));
	//		std::copy(other.coef[i].begin(), other.coef[i].end(), this->coef[i].begin());

	//		std::cout << this->coef[i] << std::endl;
	//	}

	//	std::cout << "exit copy assignment" << std::endl;

	//	return *this;
	//}

	//CubicSpline(CubicSpline&& other) :
	//	n{ other.n },
	//	left_boundary{ other.left_boundary },
	//	right_boundary{ other.right_boundary },
	//	points{other.points}
	//{
	//	std::cout << "in move constructor" << std::endl;
	//	coef.reserve(other.coef.size());

	//	for (auto i = 0; i < other.coef.size(); i++) {
	//		//std::unique_ptr<std::array<double, 4>> task = std::make_unique<std::array<double, 4>>(arr);
	//		//coef.push_back(*std::move(task));
	//		coef[i] = std::move(other.coef[i]);
	//	}

	//	other.coef = std::vector<std::array<double, 4>>();
	//	other.points = std::vector<double>();
	//	other.n = 0;
	//	other.left_boundary = 0;
	//	other.right_boundary = 0;

	//	std::cout << "exit move constructor" << std::endl;
	//}

	//CubicSpline(const CubicSpline& other) :
	//	n{ other.n },
	//	left_boundary{ other.left_boundary },
	//	right_boundary{ other.right_boundary },
	//	points{other.points}
	//{
	//	std::cout << "in copy constructor" << std::endl;
	//	coef.reserve(other.coef.size());

	//	for (auto i = 0; i < other.coef.size(); i++) {
	//		//std::unique_ptr<std::array<double, 4>> task = std::make_unique<std::array<double, 4>>(arr);
	//		//coef.push_back(*std::move(task));
	//		std::copy(other.coef[i].begin(), other.coef[i].end(), coef[i].begin());

	//		std::cout << coef[i] << std::endl;
	//	}

	//	std::cout << "exit copy constructor" << std::endl;
	//}


private:
	std::vector<double> points;
public:
	double left_boundary, right_boundary;
	int n;
	std::vector<std::array<double, 4>> coef;
};
