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
	Matrix(size_t rows_, size_t cols_);
	size_t get_rows() const;
	size_t get_cols() const;
	std::vector<T>& operator[](size_t row);
	const std::vector<T>& operator[](size_t row) const;

	Matrix<T> operator+(const Matrix<T>& other) const;
	Matrix<T> operator-(const Matrix<T>& other) const;
	Matrix<T> operator*(const Matrix<T>& other) const;
	Matrix<T> operator*(const T& scalar) const;
	Matrix<T> transpose() const;
};

template<typename T>
Matrix<T>::Matrix(size_t rows_, size_t cols_) : rows(rows_), cols(cols_)
{
	if (rows_ < 0 || cols_ < 0)
		throw std::invalid_argument("Matrix dimensions cannot be negative");
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

#endif // MATRIX_H
