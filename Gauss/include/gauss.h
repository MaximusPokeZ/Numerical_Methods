

#ifndef NUMERICAL_METHODS_GAUSS_H
#define NUMERICAL_METHODS_GAUSS_H

#include "Matrix.h"

template <typename T>
void method_gauss (Matrix<T> A, std::vector<T>& b, T& det)
{
	size_t n = A.get_rows();
	size_t count_permut = 0;

	for (size_t i = 0; i < n; ++i)
	{
		size_t max_row = i;
		for (size_t k = i + 1; k < n; ++k)
		{
			if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
			{
				max_row = k;
			}
		}

		if (i != max_row) // обмениваем строку со строкой, содержащей max элемент
		{
			std::swap(A[i], A[max_row]);
			std::swap(b[i], b[max_row]);
			++count_permut;
		}

		// прямой ход
		for (size_t j = i + 1; j < n; ++j)
		{
			if (A[j][i] != 0)
			{
				T coef = A[j][i] / A[i][i]; // нормировка делим на max_el
				for (size_t k = i; k < n; ++k)
				{
					A[j][k] -= coef * A[i][k];
				}
				b[j] -= coef * b[i];
			}
		}
	}


	// detA = (-1)^m a11 * a22 * ... * ann/ m - кол-во перестановок

	// обратный ход
	for (size_t i = n; i-- > 0;)
	{
		T sum = b[i];
		for (size_t j = i + 1; j < n; ++j)
		{
			sum -= A[i][j] * b[j];
		}
		b[i] = sum / A[i][i];
		det *= A[i][i];
	}

	det *= (count_permut % 2) ? -1 : 1;
}

template <typename T>
Matrix<T> get_reverse_matrix(Matrix<T> A)
{
	size_t n = A.get_rows();
	Matrix<T> augmented(n, 2 * n);

	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			augmented[i][j] = A[i][j];
		}
		for (size_t j = n; j < 2 * n; ++j)
		{
			augmented[i][j] = (i == (j - n)) ? 1 : 0;
		}
	}

	for (size_t i = 0; i < n; ++i)
	{
		size_t max_row = i;
		for (size_t k = i + 1; k < n; ++k)
		{
			if (std::abs(augmented[k][i]) > std::abs(augmented[max_row][i]))
			{
				max_row = k;
			}
		}

		// если ведущий элемент равен нулю, обратной матрицы не существует
		if (augmented[max_row][i] == 0) {
			throw std::runtime_error("Matrix is singular and cannot be inverted");
		}

		if (max_row != i)
		{
			for (size_t j = 0; j < 2 * n; ++j)
			{
				std::swap(augmented[i][j], augmented[max_row][j]);
			}
		}

		for (size_t k = i + 1; k < n; ++k)
		{
			T factor = augmented[k][i] / augmented[i][i];
			for (size_t j = i; j < 2 * n; ++j)
			{
				augmented[k][j] -= factor * augmented[i][j];
			}
		}
	}

	// Приведение к единичной матрице
	for (size_t i = n; i-- > 0;) {
		T divisor = augmented[i][i];
		for (size_t j = 0; j < 2 * n; ++j)
		{
			augmented[i][j] /= divisor;
		}
		for (size_t k = 0; k < i; ++k)
		{
			T factor = augmented[k][i];
			for (size_t j = 0; j < 2 * n; ++j)
			{
				augmented[k][j] -= factor * augmented[i][j];
			}
		}
	}

	Matrix<T> inverse(n, n);
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			inverse[i][j] = augmented[i][j + n]; // правая часть является обратной матрицей
		}
	}

	return inverse;
}



#endif //NUMERICAL_METHODS_GAUSS_H
