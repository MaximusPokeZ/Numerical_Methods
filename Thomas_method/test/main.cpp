#include <method.h>
#include <iomanip>

int main()
{

	try {
		std::cout << std::setprecision(6);
		size_t n;
		std::cout << "Введите размерность системы уравнений: ";
		std::cin >> n;

		std::vector<double> b(n), c(n - 1), a(n - 1), d(n);

		std::cout << "Введите коэффициенты главной диагонали a: ";
		for (size_t i = 0; i < n; ++i)
			std::cin >> b[i];

		std::cout << "Введите коэффициенты верхней диагонали b: ";
		for (size_t i = 0; i < n - 1; ++i)
			std::cin >> c[i];

		std::cout << "Введите коэффициенты нижней диагонали c: ";
		for (size_t i = 0; i < n - 1; ++i)
			std::cin >> a[i];

		std::cout << "Введите вектор правой части d: ";
		for (size_t i = 0; i < n; ++i)
			std::cin >> d[i];

		// Решение системы методом Томпсона
		double det = b[0];
		std::vector<double> result = ThomasMethod<double>::solve(b, c, a, d, det);

		// Вывод результата
		std::cout << "Решение системы:\n";
		size_t i = 0;
		for (double x : result)
		{
			std::cout << "x" << ++i << " = " << x << "\n";
		}
		std::cout << "\n";

		std::cout << "Определитель матрицы det = " << det << "\n";
	}
	catch (const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return 0;




}