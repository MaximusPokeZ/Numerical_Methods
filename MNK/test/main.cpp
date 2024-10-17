#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "Matrix.h"
#include "gauss.h"


template <typename T>
std::vector<T> solve_normal_system_mnk(const std::vector<T>& x, const std::vector<T>& y, int degree)
{
	int n = x.size();
	int m = degree + 1; // кол-во коэфф
	Matrix<T> A(m, m);
	std::vector<T> B(m, 0.0);

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			A[i][j] = 0;
			for (int k = 0; k < n; ++k)
			{
				A[i][j] += std::pow(x[k], i + j);
			}
		}
		for (int k = 0; k < n; ++k)
		{
			B[i] += std::pow(x[k], i) * y[k];
		}
	}

	// std::cout << A << "\n";

	T det = 1.0;
	method_gauss(A, B, det);

	return B;
}


template <typename T>
T evaluate_polynomial(const std::vector<T>& coeffs, T x)
{
	T result = 0;
	T power = 1;
	for (const auto& coeff : coeffs)
	{
		result += coeff * power;
		power *= x;
	}
	return result;
}


template <typename T>
T calculate_error_sum(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& coeffs)
{
	T error_sum = 0;
	for (size_t j = 0; j < x.size(); ++j)
	 {
		T F = evaluate_polynomial(coeffs, x[j]);
		error_sum += std::pow(F - y[j], 2);
	}
	return error_sum;
}

int main()
{
//	std::vector<long double> x = {0.0, 1.7, 3.4, 5.1, 6.8, 8.5};
//	std::vector<long double> y = {0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155};
	std::vector<long double> x = {-0.7, -0.4, -0.1, 0.2, 0.5, 0.8};
	std::vector<long double> y = {-0.7754, -0.41152, -0.10017, 0.20136, 0.5236, 0.9273};

	// 308 стр
	auto coeffs = solve_normal_system_mnk(x, y, 1); // F1_x = a_0 + a_1x
	std::cout << std::fixed << std::setprecision(5) << "Coefficients of the first degree polynomial: \n";
	for (size_t i = 0; i < coeffs.size(); ++i)
	{
		std::cout << "a" << i << " = " << coeffs[i] << "\n";
	}
	std::cout << "\n";

	for (size_t i = 0; i < 6; ++i)
	{
		std::cout << "F (x" << i << ") = " << evaluate_polynomial(coeffs, x[i]) << "\n";
	}
	std::cout << "\n";

	long double error_sum1 = calculate_error_sum(x, y, coeffs);
	std::cout << "Sum of squared errors for the first degree Ф: " << error_sum1 << "\n";


	coeffs = solve_normal_system_mnk(x, y, 2); // F2_x = a_0 + a_1x + a_2x^2
	std::cout << "Coefficients of the second degree polynomial: \n";

	for (size_t i = 0; i < coeffs.size(); ++i)
	{
		std::cout << "a" << i << " = " << coeffs[i] << "\n";
	}
	std::cout << "\n";

	for (size_t i = 0; i < 6; ++i)
	{
		std::cout << "F (x" << i << ") = " << evaluate_polynomial(coeffs, x[i]) << "\n";
	}
	std::cout << "\n";


	long double error_sum2 = calculate_error_sum(x, y, coeffs);
	std::cout << "Sum of squared errors for the second degree: " << error_sum2 << "\n";

	return 0;
}

0.00553 + 1.10670*x

-0.00699 + 1.10189 * x + 0.04815 * x * x