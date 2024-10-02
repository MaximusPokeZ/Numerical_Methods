#ifndef NUMERICAL_METHODS_FUNCTIONS_H
#define NUMERICAL_METHODS_FUNCTIONS_H

// http://www.dep805.ru/education/tocmsu/L10TO7.pdf

#include "Matrix.h"

template <typename T>
bool convergence_condition (const Matrix<T>& m, const long double& eps)
{
	size_t n = m.get_rows();

	for (size_t i = 0; i < n; ++i)
	{
		T sum{};
		for (size_t j = 0; j < n; ++j)
		{
			if (i != j)
				sum += std::abs(m[i][j]);
		}
		if (std::abs(m[i][i]) <= sum + eps)
			return false;
	}
	return true;
}

template <typename T>
T vector_norma(const std::vector<T>& v1, const std::vector<T>& v2)
{
	T sum = 0;
	for (size_t i = 0; i < v1.size(); ++i)
	{
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	return std::sqrt(sum);
}

template <typename T>
std::vector<T> method_simple_iterations(const Matrix<T>& m, const std::vector<T>& b, const long double& eps, size_t& count)
{
	size_t n = m.get_rows();
	std::vector<T> x(n, T{}), prev_x(n, T{});
	T norm;

	if (!convergence_condition(m, eps))
	{
		throw std::invalid_argument("The matrix does not satisfy the convergence condition!");
	}

	do
	{
		prev_x = x;
		for (size_t i = 0; i < n; ++i)
		{
			T sum = b[i];
			for (size_t j = 0; j < n; ++j)
			{
				if (i != j)
				{
					sum -= m[i][j] * prev_x[j];
				}
			}
			x[i] = sum / m[i][i];
		}
		++count;
		norm = vector_norma(x, prev_x);
	} while (norm > eps);

	return x;
}


template <typename T>
std::vector<T> method_seidel(const Matrix<T>& m, const std::vector<T>& b, const long double& eps, size_t& count)
{
	size_t n = m.get_rows();
	std::vector<T> x(n, T{}), prev_x(n, T{});
	T norm;

	if (!convergence_condition(m, eps))
	{
		throw std::invalid_argument("The matrix does not satisfy the convergence condition!");
	}

	do
	{
		prev_x = x;
		for (size_t i = 0; i < n; ++i)
		{
			T sum = T{};
			for (size_t j = 0; j < n; ++j)
			{
				if (i != j)
				{
					sum += m[i][j] * x[j];
				}
				x[i] = (b[i] - sum) / m[i][i];
			}
		}
		++count;
		norm = vector_norma(x, prev_x);
	}  while (norm > eps);

	return x;
}



#endif //NUMERICAL_METHODS_FUNCTIONS_H
