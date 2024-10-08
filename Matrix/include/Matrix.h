#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>

template <typename T>
class Matrix
{
private:
	std::vector<std::vector<T>> data;
	size_t rows, cols;

public:
	Matrix(size_t n);
	Matrix(size_t rows_, size_t cols_);
	size_t get_rows() const;
	size_t get_cols() const;
	std::vector<T>& operator[](size_t row);
	const std::vector<T>& operator[](size_t row) const;
	bool is_symmetric (T epsilon = std::numeric_limits<T>::epsilon()) const;

	Matrix<T> operator+(const Matrix<T>& other) const;
	Matrix<T> operator-(const Matrix<T>& other) const;
	Matrix<T> operator*(const Matrix<T>& other) const;
	Matrix<T> operator*(const T& scalar) const;
	Matrix<T> transpose() const;

	T frobenius_sqrt_norm() const;
	T row_norm() const;
	T col_norm() const;

};

template<typename T>
Matrix<T>::Matrix(size_t n) : rows(n), cols(n)
{
	data.resize(n, std::vector<T>(n, T{}));
	for (size_t i = 0; i < n; ++i)
	{
		data[i][i] = 1;
	}
}

template<typename T>
Matrix<T>::Matrix(size_t rows_, size_t cols_) : rows(rows_), cols(cols_)
{
	data.resize(rows_, std::vector<T>(cols_, T{}));
}

template<typename T>
size_t Matrix<T>::get_rows() const
{
	return rows;
}

template<typename T>
size_t Matrix<T>::get_cols() const
{
	return cols;
}

template<typename T>
std::vector<T>& Matrix<T>::operator[](size_t row)
{
	if (row >= rows)
		throw std::out_of_range("Index out of range");
	return data[row];
}

template<typename T>
const std::vector<T>& Matrix<T>::operator[](size_t row) const
{
	if (row >= rows)
		throw std::out_of_range("Index out of range");
	return data[row];
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
{
	if (rows != other.rows || cols != other.cols) throw std::invalid_argument("Matrix dimensions must match");
	Matrix<T> result(rows, cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			result[i][j] = data[i][j] + other[i][j];
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
{
	if (rows != other.rows || cols != other.cols) throw std::invalid_argument("Matrix dimensions must match");
	Matrix<T> result(rows, cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			result[i][j] = data[i][j] - other[i][j];
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const
{
	if (cols != other.rows) throw std::invalid_argument("Matrix multiplication dimensions mismatch");
	Matrix<T> result(rows, other.cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < other.cols; ++j)
			for (size_t k = 0; k < cols; ++k)
				result[i][j] += data[i][k] * other[k][j];
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const
{
	Matrix<T> result(rows, cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			result[i][j] = data[i][j] * scalar;
	return result;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> result(cols, rows);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			result[j][i] = data[i][j];
	return result;
}

template<typename T>
bool Matrix<T>::is_symmetric(T epsilon) const
{
	if (rows != cols) return false;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = i + 1; j < cols; ++j)
		{
			if (std::abs(data[i][j] - data[j][i]) > epsilon) return false;
		}
	}
	return true;
}


template <typename U>
std::istream& operator>>(std::istream& in, Matrix<U>& m)
{
	size_t rows = m.get_rows();
	size_t cols = m.get_cols();
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			in >> m[i][j];
		}
	}
	return in;
}

template <typename U>
std::ostream& operator<<(std::ostream& out, const Matrix<U>& m)
{
	size_t rows = m.get_rows();
	size_t cols = m.get_cols();
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			out << m[i][j] << " ";
		}
		out << "\n";
	}
	return out;
}

template <typename T>
T Matrix<T>::frobenius_sqrt_norm() const
{
	T norm = 0;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			norm += data[i][j] * data[i][j];
		}
	}
	return std::sqrt(norm);
}

template <typename T>
T Matrix<T>::row_norm() const
{
	T max_sum = 0;
	for (size_t i = 0; i < rows; ++i)
	{
		T row_sum = 0;
		for (size_t j = 0; j < cols; ++j)
		{
			row_sum += std::abs(data[i][j]);
		}
		max_sum = std::max(max_sum, row_sum);
	}
	return max_sum;
}

template <typename T>
T Matrix<T>::col_norm() const
{
	T max_sum = 0;
	for (size_t j = 0; j < cols; ++j)
	{
		T col_sum = 0;
		for (size_t i = 0; i < rows; ++i)
		{
			col_sum += std::abs(data[i][j]);
		}
		max_sum = std::max(max_sum, col_sum);
	}
	return max_sum;
}

#endif // MATRIX_H
