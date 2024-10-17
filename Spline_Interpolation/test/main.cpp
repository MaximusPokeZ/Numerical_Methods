#include <iostream>
#include <vector>
#include <iomanip>
#include "Matrix.h"
#include "gauss.h"

struct Spline
{
	long double a, b, c, d, x1, x2;
};


std::vector<Spline> build_spline(const std::vector<long double>& x, const std::vector<long double>& f)
{
	int n = x.size() - 1;

	std::vector<long double> h(n);
	for (int i = 0; i < n; ++i)
	{
		h[i] = x[i + 1] - x[i];
	}

	Matrix<long double> A(n - 1, n - 1);
	std::vector<long double> B(n - 1);

	for (int i = 1; i < n; ++i)
	{
		A[i - 1][i - 1] = 2 * (h[i - 1] + h[i]);
		if (i > 1) A[i - 1][i - 2] = h[i - 1];
		if (i < n - 1) A[i - 1][i] = h[i];
		B[i - 1] = 3 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
	}

	//std::cout << A << std::endl;
//	for (const auto & elem : B)
//	{
//		std::cout << elem << std::endl;
//	}

	long double det = 1.0;
	std::vector<long double> С(n + 1, 0.0);

	method_gauss(A, B, det);

	for (int i = 1; i < n; ++i)
	{
		С[i] = B[i - 1];
	}

	std::vector<Spline> splines(n);
	for (int i = 0; i < n; ++i)
	{
		splines[i].a = f[i];
		splines[i].b = (f[i + 1] - f[i]) / h[i] - h[i] * (С[i + 1] + 2 * С[i]) / 3;
		splines[i].c = С[i];
		splines[i].d = (С[i + 1] - С[i]) / (3 * h[i]);
		splines[i].x1 = x[i];
		splines[i].x2 = x[i + 1];
		std::cout << splines[i].a << " " << splines[i].b << " " <<splines[i].c << " " << splines[i].d << " " << splines[i].x1 << " " << splines[i].x2 <<" \n";
	}

	return splines;
}

std::pair<size_t, long double> evaluate_spline(const std::vector<Spline>& splines, long double x)
{
	size_t i = 0;
	for (const auto& spline : splines)
	{
		if (x >= spline.x1 && x <= spline.x2)
		{
			long double dx = x - spline.x1;
			return std::make_pair(i, spline.a + spline.b * dx + spline.c * dx * dx + spline.d * dx * dx * dx);
		}
		++i;
	}
	return std::make_pair(0, 0);
}

int main()
{
	std::vector<long double> x = {-0.4, -0.1, 0.2, 0.5, 0.8};
	std::vector<long double> f = {-0.41152, -0.10017, 0.20136, 0.52360, 0.92730};

//	std::vector<long double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
//	std::vector<long double> f = {0.0, 1.8415, 2.9093, 3.1411, 3.2432};

	std::vector<Spline> splines = build_spline(x, f);

	long double x_eval = 0.1;
	auto result = evaluate_spline(splines, x_eval);
	auto i = result.first;

	std::cout << std::fixed << std::setprecision(6) << "\n";
	std::cout << "Section : " << "[ " << splines[i].x1 << ", " << splines[i].x2 << "]\n";
	std::cout << "a = " << splines[i].a << "\n";
	std::cout << "b = " << splines[i].b << "\n";
	std::cout << "c = " << splines[i].c << "\n";
	std::cout << "d = " << splines[i].d << "\n";
	std::cout << "f(X*) = " << result.second << "\n";

	//выписать многочлен в общем виде

	0.5236+1.21123*(x- 0.5 ) + 0.672167*(x- 0.5 )^2 + -0.746852 * (x- 0.5 )^3

	return 0;
}
