#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

double f(double x, double y)
{
	return y * sin(x) - cos(x) * sin(x);
}

double exact_solution(double x)
{
	return cos(x) + exp(-cos(x)) - 1;
}

double runge_kutta_step(double x, double y, double h)
{
	double k1 = h * f(x, y); // вспомогательные величины
	double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
	double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
	double k4 = h * f(x + h, y + k3);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0; // приращение ф-ции
}

int main()
{
	double x0 = 0.0;
	double y0 = 0.367879;
	double h = 0.1;
	double x_end = 1.0;

	int n = static_cast<int>((x_end - x0) / h) + 1;

	std::vector<double> x_values(n);
	std::vector<double> y_euler(n);
	std::vector<double> y_runge_kutta(n);
	std::vector<double> delta_y_euler(n);
	std::vector<double> delta_y_runge_kutta(n);
	std::vector<double> errors_euler(n);
	std::vector<double> errors_runge_kutta(n);

	x_values[0] = x0;
	y_euler[0] = y0;
	y_runge_kutta[0] = y0;

	for (int i = 1; i < n; ++i)
	{
		x_values[i] = x0 + i * h;

		double y_prev = y_euler[i - 1];
		y_euler[i] = y_prev + h * f(x_values[i - 1], y_prev);
		delta_y_euler[i] = y_euler[i] - y_prev;

		y_runge_kutta[i] = runge_kutta_step(x_values[i - 1], y_runge_kutta[i - 1], h);
		delta_y_runge_kutta[i] = y_runge_kutta[i] - y_runge_kutta[i - 1];
	}

	std::cout << std::fixed << std::setprecision(8);
	std::cout << " Euler's method\n";
	std::cout << " x          y_euler       Δy_k (Euler)  y_exact      ε_euler\n";
	std::cout << "-------------------------------------------------------------------\n";
	for (int i = 0; i < n; ++i)
	{
		double y_exact_value = exact_solution(x_values[i]);
		errors_euler[i] = fabs(y_exact_value - y_euler[i]);

		std::cout << x_values[i] << "  "
				  << y_euler[i] << "  "
				  << (i == 0 ? 0.0 : delta_y_euler[i]) << "  "
				  << y_exact_value << "  "
				  << errors_euler[i] << "\n";
	}

	std::cout << "\n Runge-kutta method\n";
	std::cout << " x          y_runge_kutta  Δy_k (Runge)  y_exact      ε_runge_kutta\n";
	std::cout << "-------------------------------------------------------------------\n";
	for (int i = 0; i < n; ++i)
	{
		double y_exact_value = exact_solution(x_values[i]);
		errors_runge_kutta[i] = fabs(y_exact_value - y_runge_kutta[i]);

		std::cout << x_values[i] << "  "
				  << y_runge_kutta[i] << "  "
				  << (i == 0 ? 0.0 : delta_y_runge_kutta[i]) << "  "
				  << y_exact_value << "  "
				  << errors_runge_kutta[i] << "\n";
	}

	std::cout << std::fixed << std::setprecision(5);
	std::cout << "\n Solution by Euler's method\n";
	std::cout << "-----------------------------\nxk: ";
	for (int i = 0; i < n; ++i)
	{
		std::cout << x_values[i] << "  ";
	}
	std::cout << "\nyk: ";
	for (int i = 0; i < n; ++i)
	{
		std::cout << y_euler[i] << "  ";
	}

	std::cout << "\n\n Solution by Runge-Kutta method\n";
	std::cout << "-----------------------------\nxk: ";
	for (int i = 0; i < n; ++i)
	{
		std::cout << x_values[i] << "  ";
	}
	std::cout << "\nyk: ";
	for (int i = 0; i < n; ++i)
	{
		std::cout << y_runge_kutta[i] << "  ";
	}

	return 0;
}
