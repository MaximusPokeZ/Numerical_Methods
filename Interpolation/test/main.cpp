#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>


// четвертая производная arcsin(x)
double fourth_derivative_arcsin(double x)
{
	return -(6 * std::pow(x, 3) + 9 * x) / (std::sqrt(1 - std::pow(x, 2)) * (std::pow(x, 6) - 3 * std::pow(x, 4) + 3 * std::pow(x, 2) - 1));
}


double find_max_fourth_derivative(const double & a, const double & b)
{
	double max_value = 0.0;
	double x = a;
	for (; x <= b; x += 0.05)
	{
		double value = std::abs(fourth_derivative_arcsin(x));
		if (value > max_value)
		{
			max_value = value;
		}
	}
	return max_value;
}


double divided_difference(const std::vector<double>& x, const std::vector<double>& y, int i, int j)
{
	if (i == j)
	{
		return y[i];
	}
	else
	{
		return (divided_difference(x, y, i + 1, j) - divided_difference(x, y, i, j - 1)) / (x[j] - x[i]);
	}
}

double newton_polynomial(const std::vector<double>& x, const std::vector<double>& y, double x_star)
{
	double result = y[0];
	double product = 1.0;

	for (size_t i = 1; i < x.size(); ++i)
	{
		product *= (x_star - x[i - 1]);
		result += divided_difference(x, y, 0, i) * product;
	}


	return result;
}

double lagrange_polynomial(const std::vector<double>& x, const std::vector<double>& y, double x_star)
{
	double result = 0.0;

	for (size_t i = 0; i < x.size(); ++i)
	{
		double term = y[i];
		for (size_t j = 0; j < x.size(); ++j)
		{
			if (j != i)
			{
				term *= (x_star - x[j]) / (x[i] - x[j]);
			}
		}
		result += term;
	}

	return result;
}

double error_estimate(const std::vector<double>& x, double max_derivative, double x_star)
{
	double omega = 1.0;
	for (double xi : x)
	{
		omega *= (x_star - xi);
	}
	return std::abs(max_derivative * omega / 24); // (n+1)! = 4!, т.к  3го порядка многочлен
}

double omega_derivative(const std::vector<double>& xi, int i)
{
	double result = 1.0;
	for (size_t j = 0; j < xi.size(); ++j)
	{
		if (j != i)
		{
			result *= (xi[i] - xi[j]);
		}
	}
	return result;
}

std::string lagrange_dop_(const std::vector<double>& xi, int i)
{
	std::string basis = "";
	for (size_t j = 0; j < xi.size(); ++j)
	{
		if (j != i)
		{
			basis += "(x - " + std::to_string(xi[j]) + ")";
		}
	}
	return basis;
}


std::string build_lagrange_polynomial(const std::vector<double>& xi, const std::vector<double>& yi)
{
	std::string polynomial = "";
	for (size_t i = 0; i < xi.size(); ++i)
	{
		double omega_d4 = omega_derivative(xi, i);
		double coeff = yi[i] / omega_d4;
		std::string term = std::to_string(coeff) + " * " + lagrange_dop_(xi, i);
		if (i > 0 && coeff >= 0)
		{
			polynomial += " + " + term;
		}
		else
		{
			polynomial += " " + term;
		}
	}
	return polynomial;
}

std::string build_newton_polynomial(const std::vector<double>& xi, const std::vector<double>& yi)
{
	std::ostringstream polynomial;
	polynomial << std::fixed << std::setprecision(2);

	std::vector<double> coefficients;
	for (size_t i = 0; i < yi.size(); ++i)
	{
		coefficients.push_back(divided_difference(yi, xi, 0, i));
	}

	polynomial << coefficients[0];

	for (size_t i = 1; i < coefficients.size(); ++i)
	{
		std::string term = std::to_string(coefficients[i]);
		for (size_t j = 0; j < i; ++j)
		{
			term += "(x - " + std::to_string(xi[j]) + ")";
		}
		polynomial << " + " << term;
	}

	return polynomial.str();
}

int main()
{
	double a = -0.4;
	double b = 0.5;
	std::vector<double> x = {-0.4, -0.1, 0.2, 0.5};
	//std::vector<double> x = {-0.4, 0, 0.2, 0.5};
	std::vector<double> y;
	for (double xi : x)
	{
		y.push_back(std::asin(xi)); // y = arcsin(x)
	}

	std::string lagrange_polynomial_str = build_lagrange_polynomial(x, y);
	std::cout << "L_3(x) = " << lagrange_polynomial_str << std::endl;

	std::string newton_polynomial_str = build_newton_polynomial(x, y);
	std::cout << "P_3(x) = " << newton_polynomial_str << std::endl;

	double x_star = 0.1;

	double lagrange_value = lagrange_polynomial(x, y, x_star);
	double newton_value = newton_polynomial(x, y, x_star);

	double actual_value = std::asin(x_star);

	double lagrange_error = std::abs(actual_value - lagrange_value);
	double newton_error = std::abs(actual_value - newton_value);

	// Оценка погрешности по формуле (3.9)
	double max_derivative = find_max_fourth_derivative(a, b);
	double error_bound = error_estimate(x, max_derivative, x_star);

	std::cout << std::fixed << std::setprecision(5);

	std::cout << "Lagrange interpolation value at x* = " << x_star << ": ---> " << lagrange_value << std::endl;
	std::cout << "Newton interpolation value at x* = " << x_star << ": ---> " << newton_value << std::endl;
	std::cout << "Actual value of arcsin(" << x_star << "): " << actual_value << std::endl;

	std::cout << "Lagrange interpolation error: " << lagrange_error << std::endl;
	std::cout << "Newton interpolation error: " << newton_error << std::endl;
	std::cout << "Estimated error bound (3.9): " << error_bound << std::endl;

	return 0;
}
