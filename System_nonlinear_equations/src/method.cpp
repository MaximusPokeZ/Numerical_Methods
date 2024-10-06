#include "../include/method.h"

// first eq
long double f1 (long double x1, long double x2)
{
	return x1 * x1 + x2 * x2 - 9;
}

long double dx1_f1(long double x1, long double x2)
{
	return 2 * x1;
}

long double dx2_f1(long double x1, long double x2)
{
	return 2 * x2;
}

// second eq
long double f2 (long double x1, long double x2)
{
	return x1 - std::exp(x2) + 3;
}

long double dx1_f2(long double x1, long double x2)
{
	return 1;
}

long double dx2_f2(long double x1, long double x2)
{
	return -std::exp(x2);
}

long double determinant_2x2(const Matrix<long double>& J)
{
	return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

size_t newton_method_system(const long double eps, long double& x1, long double& x2,
						  Matrix<std::function<long double(long double, long double)>> J,
						  Matrix<std::function<long double(long double, long double)>> A1,
						  Matrix<std::function<long double(long double, long double)>> A2)
{

	long double prev_x1, prev_x2;
	size_t iteration = 0;
	Matrix<long double> J_val(2, 2);
	Matrix<long double> A1_val(2, 2);
	Matrix<long double> A2_val(2, 2);

	do
	{
		prev_x1 = x1;
		prev_x2 = x2;

		J_val[0][0] = J[0][0](prev_x1, prev_x2); J_val[0][1] = J[0][1](prev_x1, prev_x2);
		J_val[1][0] = J[1][0](prev_x1, prev_x2); J_val[1][1] = J[1][1](prev_x1, prev_x2);

		A1_val[0][0] = A1[0][0](prev_x1, prev_x2); A1_val[0][1] = A1[0][1](prev_x1, prev_x2);
		A1_val[1][0] = A1[1][0](prev_x1, prev_x2); A1_val[1][1] = A1[1][1](prev_x1, prev_x2);

		A2_val[0][0] = A2[0][0](prev_x1, prev_x2); A2_val[0][1] = A2[0][1](prev_x1, prev_x2);
		A2_val[1][0] = A2[1][0](prev_x1, prev_x2); A2_val[1][1] = A2[1][1](prev_x1, prev_x2);

		// Определитель Якобиана
		long double det = determinant_2x2(J_val);
		if (std::abs(det) < eps)
		{
			throw std::runtime_error("Якобиан вырожден!");
		}

		x1 = prev_x1 - determinant_2x2(A1_val) / det;
		x2 = prev_x2 - determinant_2x2(A2_val) / det;

		iteration++;

	} while (std::abs(x1 - prev_x1) > eps || std::abs(x2 - prev_x2) > eps);

	return iteration;
}


long double phi1 (long double x1, long double x2)
{
	return std::sqrt(9 - x2 * x2);
}

long double phi1_dx1 (long double x1, long double x2)
{
	return 0;
}

long double phi1_dx2 (long double x1, long double x2)
{
	return -x2 / std::sqrt(9 - x2 * x2);
}

long double phi1_m (long double x1, long double x2)
{
	return -std::sqrt(9 - x2 * x2);
}

long double phi1_dx1_m (long double x1, long double x2)
{
	return 0;
}

long double phi1_dx2_m (long double x1, long double x2)
{
	return x2 / std::sqrt(9 - x2 * x2);
}

long double phi2 (long double x1, long double x2)
{
	return log(x1 + 3);
}

long double phi2_dx1 (long double x1, long double x2)
{
	return 1 / (x1 + 3);
}

long double phi2_dx2 (long double x1, long double x2)
{
	return 0;
}

size_t simple_iteration_method(long double eps, long double& x1, long double& x2,
							 std::function<long double(long double, long double)> phi1,
							 std::function<long double(long double, long double)> phi2, const long double& q)
 {
	long double prev_x1, prev_x2;
	size_t iteration = 0;

	do
	{
		prev_x1 = x1;
		prev_x2 = x2;

		x1 = phi1(prev_x1, prev_x2);
		x2 = phi2(prev_x1, prev_x2);

		iteration++;

	} while ((q / (1 - q)) * std::max(std::abs(x1 - prev_x1), std::abs(x2 - prev_x2)) > eps);

	 return iteration;
}







