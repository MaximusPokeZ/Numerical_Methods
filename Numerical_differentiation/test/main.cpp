#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>


size_t find_index(long double eps, long double x_star, const std::vector<long double>& vector_x)
{
	size_t n = vector_x.size();
	for (size_t i = 0; i < n; ++i)
	{
		if (std::abs(x_star - vector_x[i]) < eps)
		{
			return i;
		}
	}
	throw std::logic_error("x_star not found in vector x");
}


long double leftist_diff(long double eps, size_t ind, const std::vector<long double>& vector_y, const std::vector<long double>& vector_x)
{
	if (ind == 0)
		throw std::logic_error("The left derivative cannot be calculated for the first element");

	long double dx = vector_x[ind] - vector_x[ind - 1];
	if (std::abs(dx) < eps)
		throw std::runtime_error("Division by zero when calculating left diff!");

	return (vector_y[ind] - vector_y[ind - 1]) / dx;
}


long double rightist_diff(long double eps, size_t ind, const std::vector<long double>& vector_y, const std::vector<long double>& vector_x)
{
	if (ind >= vector_x.size() - 1)
		throw std::logic_error("The right derivative cannot be calculated for the last element");

	long double dx = vector_x[ind + 1] - vector_x[ind];
	if (std::abs(dx) < eps)
		throw std::runtime_error("Division by zero when calculating right diff!");

	return (vector_y[ind + 1] - vector_y[ind]) / dx;
}


long double centre_diff(long double eps, size_t ind, long double x_star, const std::vector<long double>& vector_y, const std::vector<long double>& vector_x)
{
	if (ind >= vector_x.size() - 2)
		throw std::logic_error("The central derivative cannot be calculated for the last two elements");

	long double dx1 = vector_x[ind + 1] - vector_x[ind];
	long double dx2 = vector_x[ind + 2] - vector_x[ind + 1];
	long double dx = vector_x[ind + 2] - vector_x[ind];

	if (std::abs(dx1) < eps || std::abs(dx2) < eps || std::abs(dx) < eps)
		throw std::runtime_error("Division by zero when calculating center diff!");

	long double term1 = (vector_y[ind + 1] - vector_y[ind]) / dx1;
	long double term2 = (vector_y[ind + 2] - vector_y[ind + 1]) / dx2;

	return term1 + (term2 - term1) / dx * (2 * x_star - vector_x[ind] - vector_x[ind + 1]);
}


long double second_diff(long double eps, size_t ind, long double x_star, const std::vector<long double>& vector_y, const std::vector<long double>& vector_x)
{
	if (ind >= vector_x.size() - 2)
		throw std::logic_error("The second derivative cannot be calculated for the last two elements");

	long double dx1 = vector_x[ind + 1] - vector_x[ind];
	long double dx2 = vector_x[ind + 2] - vector_x[ind + 1];
	long double dx = vector_x[ind + 2] - vector_x[ind];

	if (std::abs(dx1) < eps || std::abs(dx2) < eps || std::abs(dx) < eps)
		throw std::runtime_error("Division by zero when calculating second diff!");

	long double term1 = (vector_y[ind + 1] - vector_y[ind]) / dx1;
	long double term2 = (vector_y[ind + 2] - vector_y[ind + 1]) / dx2;

	return 2 * (term2 - term1) / dx;
}

int main()
{
	std::vector<long double> vector_x = {-1.0, 0.0, 1.0, 2.0, 3.0};
	std::vector<long double> vector_y = {-0.7854, 0.0, 0.7854, 1.1071, 1.249};
	long double eps = 0.00001;
	long double x_star = 1.0;

	try
	{
		size_t index = find_index(eps, x_star, vector_x);
		long double left = 0, right = 0, centre = 0, second = 0;

		if (index > 0)
			left = leftist_diff(eps, index, vector_y, vector_x);

		if (index < vector_x.size() - 1)
			right = rightist_diff(eps, index, vector_y, vector_x);

		if (index > 0 && index < vector_x.size() - 2)
		{
			centre = centre_diff(eps, index - 1, x_star, vector_y, vector_x);
			second = second_diff(eps, index - 1, x_star, vector_y, vector_x);
		}

		std::cout << std::fixed << std::setprecision(8);
		if (index > 0) std::cout << "leftist differential: " << left << "\n";
		if (index < vector_x.size() - 1) std::cout << "rightist differential: " << right << "\n";
		if (index > 0 && index < vector_x.size() - 2) std::cout << "centre differential: " << centre << "\n";
		if (index > 0 && index < vector_x.size() - 2) std::cout << "second differential: " << second << "\n";
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
	}

	return 0;
}
