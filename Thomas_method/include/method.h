#ifndef NUMERICAL_METHODS_METHOD_H
#define NUMERICAL_METHODS_METHOD_H

#include <Matrix.h>

template <typename T>
class ThomasMethod
{
public:
	static std::vector<T> solve(const std::vector<T>& b,
								const std::vector<T>& c,
								const std::vector<T>& a,
								const std::vector<T>& d,
								T& det
								);
};

template <typename T>
std::vector<T> ThomasMethod<T>::solve(const std::vector<T>& b,
									  const std::vector<T>& c,
									  const std::vector<T>& a,
									  const std::vector<T>& d,
									  T& det
									  )
{
	size_t n = b.size();

	if (c.size() != n - 1 || a.size() != n - 1 || d.size() != n)
		throw std::invalid_argument("Неправильные размеры векторов для трехдиагональной матрицы");

	for (size_t i = 0; i < n; ++i)
	{
		T up_down_diag_sum = 0;
		if (i < n - 1) up_down_diag_sum += std::abs(c[i]);
		if (i > 0) up_down_diag_sum += std::abs(a[i - 1]);

		if (std::abs(b[i]) < up_down_diag_sum) {
			throw std::invalid_argument("Матрица не удовлетворяет условию диагонального преобладания");
		}
	}

	std::vector<T> P(n), Q(n);

	// Прямой ход
	P[0] = -c[0] / b[0];
	Q[0] = d[0] / b[0];
	for (size_t i = 1; i < n; ++i)
	{
		T del = b[i] + a[i - 1] * P[i - 1];
		if (del == 0) throw std::runtime_error("Деление на ноль при прямом ходе в алгоритме Томаса");

		P[i] = (i < n - 1) ? -c[i] / del : 0;  // На последнем шаге P(n-1) = 0
		Q[i] = (d[i] - a[i - 1] * Q[i - 1]) / del;

		det *= del;

		if (std::abs(P[i]) > 1) throw std::invalid_argument("Условие |Pi| <= 1 нарушено");
	}

	// Обратный ход
	std::vector<T> x(n);
	x[n - 1] = Q[n - 1]; // P(n-1) = 0

	for (size_t i = n - 2; i < n; --i)
	{
		x[i] = P[i] * x[i + 1] + Q[i];
	}

	return x;
}

#endif //NUMERICAL_METHODS_METHOD_H
