#include "../include/methods.h"

long double dichotomy (long double eps, func f, long double a, long double b)
{
	if (f(a) * f(b) > 0) throw std::logic_error("The same sign at the ends of segment!");

	long double mid;
	while ((b - a) > 2 * eps)
	{
		mid = (a + b) * 0.5;

		if (std::abs(f(mid)) < eps) return mid;

		if (f(a) * f(mid) < 0) b = mid;
		else
			a = mid;
	}

	return (a + b) * 0.5;
}

long double newton_method(long double eps, long double a, long double b, std::function<long double(long double)> func,
						  std::function<long double(long double)> dfunc, std::function<long double(long double)> ddfunc)
{
	long double x, prev_x;

	if (func(a) * func(b) > 0)
	{
		throw std::logic_error("The same sign at the ends of segment!");
	}

	if (func(a) * ddfunc(a) > 0)
	{
		x = a;
	}
	else if (func(b) * ddfunc(b) > 0)
	{
		x = b;
	}
	else
	{
		x = (a + b) * 0.5;
	}

	long double df;
	do
	{

		prev_x = x;

		df = dfunc(x);

		if (std::abs(df) < eps) throw std::runtime_error("Division by zero!");

		x = x - func(x) / df;

	}
	while (std::abs(prev_x - x) > eps);

	return x;
}


long double secant_method(long double eps, long double a, long double b, std::function<long double(long double)> func, std::function<long double(long double)> ddfunc)
{
	long double x, prev_x, next_x;

	if (func(a) * func(b) > 0)
	{
		throw std::logic_error("The same sign at the ends of segment!");
	}

	if (func(a) * ddfunc(a) > 0)
	{
		prev_x = a;
	}
	else if (func(b) * ddfunc(b) > 0)
	{
		prev_x = b;
	}
	else
	{
		prev_x = (a + b) * 0.5;
	}

	x = prev_x + eps;

	long double fx;
	long double fxp;
	do
	{
		fx = func(x);
		fxp = func(prev_x);

		next_x = x - (fx * (x - prev_x)) / (fx - fxp);

		prev_x = x;
		x = next_x;
	}
	while (std::abs(prev_x - x) > eps);

	return x;
}

long double simple_iter (const std::function<long double(long double)>& phi_f, const std::function<long double(long double)>& d_phi_f, long double a, long double b, long double eps)
{
	long double q;
	long double prev_x = (a + b) * 0.5, x = prev_x;

	do
	{
		q = std::abs(d_phi_f(x));
		if (q >= 1) throw std::runtime_error("The convergence condition is not satisfied!");

		prev_x = x;
		x = phi_f(prev_x);

	} while(((q / (1 - q)) * std::abs(x - prev_x)) > eps);

	return x;
}


long double func8 (long double x)
{
	return log(x + 1.0) - 2.0 * x * x + 1.0;
}

long double dfunc8(long double x)
{
	return 1.0 / (x + 1.0) - 4.0 * x;
}

long double ddfunc8(long double x)
{
	return -1.0 / ((x + 1.0) * (x + 1.0)) - 4.0;
}

long double phi_x_p(long double x)
{
	return std::sqrt((log(x + 1) + 1) * 0.5);
}

long double dphi_x_p (long double x)
{
	return 1 / (2 * sqrt(2) * (x + 1) * sqrt(log(x + 1) + 1));
}
