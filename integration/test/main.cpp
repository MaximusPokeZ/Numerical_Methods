#include <iostream>
#include <cmath>

double func(double x)
{
	return 1.0 / (x * x + 4);
}

double second_derivative(double x)
{
	return (6 * x * x - 8) / (pow(x, 6) + 12 * pow(x, 4) + 48 * x * x + 64);
}

double fourth_derivative(double x)
{
	double numerator = 120 * pow(x, 4) - 960 * pow(x, 2) + 384;
	double denominator = pow(x, 10) + 20 * pow(x, 8) + 160 * pow(x, 6) + 640 * pow(x, 4) + 1280 * pow(x, 2) + 1024;
	return numerator / denominator;
}

double trapezoidal_method(double x0, double x1, double h)
{
	double sum = 0.0;
	for (double xi  = x0; xi <= x1; xi += h)
	{
		sum += func(xi);
	}
	sum += (func(x0) + func(x1)) * 0.5;
	sum *= h;
	return sum;
}


double rectangle_method(double x0, double x1, double h)
{
	double sum = 0.0;
	for (double xi  = x0; xi <= x1 - h; xi += h)
	{
		sum += func((xi + xi + h) * 0.5);
	}
	sum *= h;
	return sum;
}


double simpson_method(double x0, double x1, double h)
{
	double sum = 0.0;
	int n = (x1 - x0) / h;


	sum += (func(x0) + func(x1));
	for (int i = 1; i < n; ++i) {
		double xi = x0 + i * h; // - узел разностной сетки
		if (i % 2 == 0)
		{
			sum += 2 * func(xi);
		}
		else
		{
			sum += 4 * func(xi);
		}
	}
	sum *= h / 3.0;
	return sum;
}

double get_R_for_rectangle(long double h, long double x0, long double x1)
{
	return h * h * (x1 - x0) * std::abs(second_derivative(0)) / 24;
}

long double get_R_for_trapezoidal(long double h, long double x0, long double x1)
{
	return h * h * (x1 - x0) * std::abs(second_derivative(0)) / 12;
}

long double get_R_for_simpson(long double h, long double x0, long double x1)
{
	return h * h * h * h * (x1 - x0) * std::abs(fourth_derivative(0)) / 180;
}

long double runge_romberg_richardson(long double F_h, long double F_kh, long double k, long double p)
{
	return F_h + (F_h - F_kh)/(std::pow(k, p) - 1);
}


int main()
{
	long double h1 = 1.0;
	long double h2 = 0.5;
	long double x0 = -2.0;
	long double x1 = 2.0;

	double R_rect_h1 = get_R_for_rectangle(h1, x0, x1);
	double R_trap_h1 = get_R_for_trapezoidal(h1, x0, x1);
	double R_simp_h1 = get_R_for_simpson(h1, x0, x1);

	double R_rect_h2 = get_R_for_rectangle(h2, x0, x1);
	double R_trap_h2 = get_R_for_trapezoidal(h2, x0, x1);
	double R_simp_h2 = get_R_for_simpson(h2, x0, x1);

	// Выводим результаты
	std::cout << "Residual terms for h1 = 1.0:" << std::endl;
	std::cout << "Rectangle method: R(h1 = 1.0) <= " << R_rect_h1 << std::endl;
	std::cout << "Trapezoidal method: R(h1 = 1.0) <= " << R_trap_h1 << std::endl;
	std::cout << "Simpson method: R(h1 = 1.0) <= " << R_simp_h1 << std::endl;

	std::cout << "\nResidual terms for h2 = 0.5:" << std::endl;
	std::cout << "Rectangle method: R(h2 = 0.5) <= " << R_rect_h2 << std::endl;
	std::cout << "Trapezoidal method: R(h2 = 0.5) <= " << R_trap_h2 << std::endl;
	std::cout << "Simpson method: R(h2 = 0.5) <= " << R_simp_h2 << std::endl;

	long double F_rect_h1 = rectangle_method(x0, x1, h1);
	long double F_rect_h2 = rectangle_method(x0, x1, h2);

	long double F_trap_h1 =  trapezoidal_method(x0, x1, h1);
	long double F_trap_h2 =  trapezoidal_method(x0, x1, h2);

	long double F_simp_h1 = simpson_method(x0, x1, h1);
	long double F_simp_h2 = simpson_method(x0, x1, h2);

	long double RRR_rect = runge_romberg_richardson(F_rect_h1, F_rect_h2, h1 / h2, 2);
	long double RRR_trap = runge_romberg_richardson(F_trap_h1, F_trap_h2, h1 / h2, 2);
	long double RRR_simp = runge_romberg_richardson(F_simp_h1, F_simp_h2, h1 / h2, 4);

	long double abs_error_rect = std::abs(F_rect_h2 - RRR_rect);
	long double abs_error_trap = std::abs(F_trap_h2 - RRR_trap);
	long double abs_error_simp = std::abs(F_simp_h2 - RRR_simp);

	std::cout << "\nResult for h1 = 1.0:" << std::endl;
	std::cout << "Rectangle method: " << F_rect_h1 << std::endl;
	std::cout << "Trapezoidal method: " << F_trap_h1 << std::endl;
	std::cout << "Simpson method: " << F_simp_h1 << std::endl;

	std::cout << "\nResult for h2 = 0.5:" << std::endl;
	std::cout << "Rectangle method: " << F_rect_h2 << std::endl;
	std::cout << "Trapezoidal method: " << F_trap_h2 << std::endl;
	std::cout << "Simpson method: " << F_simp_h2 << std::endl;

	std::cout << "\nRunge-Romberg-Richardson corrected values:" << std::endl; // p - порядок точности k - коэффициент отношения шагов
	std::cout << "Rectangle method: " << RRR_rect << std::endl;
	std::cout << "Trapezoidal method: " << RRR_trap << std::endl;
	std::cout << "Simpson method: " << RRR_simp << std::endl;

	std::cout << "\nAbsolute errors:" << std::endl;
	std::cout << "Rectangle method: " << abs_error_rect << std::endl;
	std::cout << "Trapezoidal method: " << abs_error_trap << std::endl;
	std::cout << "Simpson method: " << abs_error_simp << std::endl;

	return 0;
}
