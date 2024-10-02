#include "functions.h"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	try {

		if (argc != 2) {
			std::cout << argc;
			std::cerr << "Incorrect count of parameters";
			return 1;
		}

		std::ifstream infile(argv[1]);

		if (!infile.is_open()) {
			std::cerr << "Incorrect file input";
			return 1;
		}

		size_t n;
		infile >> n;

		Matrix<long double> matrix(n, n);
		std::vector<long double> b(n);

		infile >> matrix;
		for (size_t i = 0; i < n; ++i) {
			infile >> b[i];
		}

		long double eps;
		infile >> eps;

		infile.close();

		std::cout << "Matrix A:\n";
		std::cout << matrix;

		std::cout << "Vector b:\n";
		for (long double elem: b) {
			std::cout << elem << "\n";
		}

		size_t count = 0;
		std::vector<long double> res = method_simple_iterations(matrix, b, eps, count);

		std::cout << "Result:\n";
		for (long double elem: res) {
			std::cout << elem << "\n";
		}

		std::cout << "\nCount of iterations " << count;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}


	return 0;
}