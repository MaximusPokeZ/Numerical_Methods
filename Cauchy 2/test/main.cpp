#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

double f(double x, double y, double z)
{
	return 4 * x * z - (4 * x * x - 3) * y + exp(x * x);
}

double exact_solution(double x)
{
	return (exp(x) + exp(-x) - 1) * exp(x * x);
}

void euler_method(std::vector<double>& x_values, std::vector<double>& y_values, std::vector<double>& z_values, double h)
{
	for (size_t i = 1; i < x_values.size(); ++i)
	{
		double x = x_values[i - 1];
		double y = y_values[i - 1];
		double z = z_values[i - 1];

		y_values[i] = y + h * z;
		z_values[i] = z + h * f(x, y, z);

		x_values[i] = x + h;
	}
}


void runge_kutta_4(std::vector<double>& x_values, std::vector<double>& y_values, std::vector<double>& z_values, double h, size_t n)
{
	for (size_t i = 1; i < n; ++i)
	{
		double x = x_values[i - 1];
		double y = y_values[i - 1];
		double z = z_values[i - 1];

		double k1_y = h * z;
		double L1_z = h * f(x, y, z);

		double k2_y = h * (z + 0.5 * L1_z);
		double L2_z = h * f(x + 0.5 * h, y + 0.5 * k1_y, z + 0.5 * L1_z);

		double k3_y = h * (z + 0.5 * L2_z);
		double L3_z = h * f(x + 0.5 * h, y + 0.5 * k2_y, z + 0.5 * L2_z);

		double k4_y = h * (z + L3_z);
		double L4_z = h * f(x + h, y + k3_y, z + L3_z);


		// значения и приращения функций
		y_values[i] = y + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
		z_values[i] = z + (L1_z + 2 * L2_z + 2 * L3_z + L4_z) / 6;

		x_values[i] = x + h;
	}
}

void adams_bashforth_4(std::vector<double>& x_values, std::vector<double>& y_values, std::vector<double>& z_values, double h)
{
	size_t n = x_values.size();

	runge_kutta_4(x_values, y_values, z_values, h, 4); 

	for (size_t i = 4; i < n; ++i)
	{
		double y0 = y_values[i - 1];
		double z0 = z_values[i - 1];

		double y_pred = y0 + h * (55 * z_values[i - 1] - 59 * z_values[i - 2] + 37 * z_values[i - 3] - 9 * z_values[i - 4]) / 24;

		double z_pred = z0 + h * (55 * f(x_values[i - 1], y_values[i - 1], z_values[i - 1])
								  - 59 * f(x_values[i - 2], y_values[i - 2], z_values[i - 2])
								  + 37 * f(x_values[i - 3], y_values[i - 3], z_values[i - 3])
								  - 9 * f(x_values[i - 4], y_values[i - 4], z_values[i - 4])) / 24;

		// корректор
		double f_pred = f(x_values[i - 1] + h, y_pred, z_pred);

		y_values[i] = y0 + h * (9 * z_pred + 19 * z_values[i - 1] - 5 * z_values[i - 2] + z_values[i - 3]) / 24;

		z_values[i] = z0 + h * (9 * f_pred + 19 * f(x_values[i - 1], y_values[i - 1], z_values[i - 1])
								- 5 * f(x_values[i - 2], y_values[i - 2], z_values[i - 2])
								+ f(x_values[i - 3], y_values[i - 3], z_values[i - 3])) / 24;

		x_values[i] = x_values[i - 1] + h;
	}
}

int main()
{
	double x0 = 0.0, x_end = 1.0;
	double y0 = 1.0, y0_d = 0.0;
	double h = 0.1;
	int n = static_cast<int>((x_end - x0) / h) + 1;

	std::vector<double> x_values(n), y_euler(n), z_euler(n), y_adams(n), z_adams(n), y_runge_kutta(n), z_runge_kutta(n);
	std::vector<double> y_exact(n), errors_euler(n), errors_adams(n), errors_runge_kutta(n);

	x_values[0] = x0;
	y_euler[0] = y_adams[0] = y_runge_kutta[0] = y0;
	z_euler[0] = z_adams[0] = z_runge_kutta[0] = y0_d;


	euler_method(x_values, y_euler, z_euler, h);

	runge_kutta_4(x_values, y_runge_kutta, z_runge_kutta, h, x_values.size());

	adams_bashforth_4(x_values, y_adams, z_adams, h);

	std::cout << "Euler's method:\n" << std::setprecision(8);
	std::cout << " x          y_numerical    y_exact      error\n";
	std::cout << "--------------------------------------------\n";
	for (int i = 0; i < n; ++i)
	{
		double y_exact_value = exact_solution(x_values[i]);
		errors_euler[i] = fabs(y_exact_value - y_euler[i]);
		std::cout << x_values[i] << "  " << y_euler[i] << "  " << y_exact_value << "  " << errors_euler[i] << "\n";
	}

	std::cout << "\nRunge-Kutta method of the 4th order:\n";
	std::cout << " x          y_numerical    y_exact      error\n";
	std::cout << "--------------------------------------------\n";
	for (int i = 0; i < n; ++i)
	{
		double y_exact_value = exact_solution(x_values[i]);
		errors_runge_kutta[i] = fabs(y_exact_value - y_runge_kutta[i]);
		std::cout << x_values[i] << "  " << y_runge_kutta[i] << "  " << y_exact_value << "  " << errors_runge_kutta[i] << "\n";
	}

	std::cout << "\nAdams-Bashforth method:\n";
	std::cout << " x          y_numerical    y_exact      error\n";
	std::cout << "--------------------------------------------\n";
	for (int i = 0; i < n; ++i)
	{
		double y_exact_value = exact_solution(x_values[i]);
		errors_adams[i] = fabs(y_exact_value - y_adams[i]);
		std::cout << x_values[i] << "  " << y_adams[i] << "  " << y_exact_value << "  " << errors_adams[i] << "\n";
	}

	return 0;
}
