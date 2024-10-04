#ifndef NUMERICAL_METHODS_ROTATE_METHOD_H
#define NUMERICAL_METHODS_ROTATE_METHOD_H

#include "Matrix.h"
#include "cmath"

// https://24calc.ru/metod-yakobi/?ysclid=m1sc1522wo827202540

template <typename T>
T not_diagonal_max_element (const Matrix<T>& m, size_t& p, size_t& q)
{
	T max_val = 0;
	size_t n = m.get_rows();

	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = i + 1; j < n; ++j)
		{
			if (std::abs(m[i][j]) > max_val)
			{
				max_val = std::abs(m[i][j]);
				p = i;
				q = j;
			}
		}
	}
	return max_val;
}


template <typename T>
void jacobi_method(Matrix<T> A, std::vector<T>& eigenvalues, Matrix<T>& J, const T& epsilon, size_t& iterations)
{
	size_t n = A.get_rows();

	if (!A.is_symmetric())
	{
		throw std::invalid_argument("The matrix does not symmetric!");
	}

	eigenvalues.resize(n);
	iterations = 0;
	size_t p = 0, q = 0;
	T max_elem = not_diagonal_max_element(A, p, q);;

	 while (max_elem > epsilon) // можно еще проверить как сумму квадратов внедиагональных элементов
	 {
		 double phi;
		 if (std::abs(A[p][p] - A[q][q]) < epsilon)
		 {
			 phi = M_PI / 4.0;
		 }
		 phi = 0.5 * std::atan(2 * A[p][q] / (A[p][p] - A[q][q]));
		 Matrix<T> P = Matrix<T>(n);
		 P[p][p] = P[q][q] = std::cos(phi);
		 P[p][q] = -std::sin(phi);
		 P[q][p] = std::sin(phi);

		 A = P.transpose() * A * P;
		 J = J * P;

		 max_elem = not_diagonal_max_element(A, p, q);;
		 ++iterations;
	 }

	 for (size_t i = 0; i < n; ++i)
	 {
		 eigenvalues[i] = A[i][i];
	 }
}

template<typename T>
std::vector<T> norma_v(const std::vector<T>& v)
{
	T norm = 0;
	size_t n = v.size();
	for (size_t i = 0; i < n; ++i)
	{
		norm += v[i] * v[i];
	}
	norm = std::sqrt(norm);

	std::vector<T> v_n = v;
	for (size_t i = 0; i < n; ++i)
	{
		v_n[i] /= norm;
	}

	return v_n;
}

// 253 стр пирумов
template <typename T>
std::pair<T, std::vector<T>> power_method(const Matrix<T>& A, const T& epsilon, size_t& iteration)
{
	size_t n = A.get_rows();
	std::vector<T> y(n, 1.0);
	std::vector<T> z(n, 0.0);
	T lambda_old = 0, lambda_new = 1;

	for (;; ++iteration)
	{
		lambda_old = (iteration) ? lambda_new : lambda_old;

		for (size_t i = 0; i < n; ++i) // z = A * y
		{
			z[i] = 0;
			for (size_t j = 0; j < n; ++j)
			{
				z[i] += A[i][j] * y[j];
			}
		}

		y = (iteration) ? norma_v(z) : y;

		lambda_new = z[0] / y[0];

		y = (iteration) ? y : norma_v(z) ;

		if (std::abs(lambda_new - lambda_old) < epsilon) break;
	}

	z = norma_v(y);

	return {lambda_new, z};
}




#endif //NUMERICAL_METHODS_ROTATE_METHOD_H
