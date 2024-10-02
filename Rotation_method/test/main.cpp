#include "rotate_method.h"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	try {

		if (argc != 2)
		{
			std::cerr << "Incorrect count of parameters";
			return 1;
		}

		std::ifstream infile(argv[1]);

		if (!infile.is_open())
		{
			std::cerr << "Incorrect file input";
			return 1;
		}

		size_t n;
		infile >> n;

		Matrix<long double> matrix(n, n);

		infile >> matrix;

		long double eps;
		infile >> eps;

		infile.close();

		std::cout << "Matrix A:\n";
		std::cout << matrix;

		Matrix<long double> J(n);
		std::vector<long double> eigenvalues(n);
		size_t iter = 0;

		jacobi_method(matrix, eigenvalues, J, eps, iter);

		std::cout << "Eigenvalues:\n";
		for (const auto& value : eigenvalues)
		{
			std::cout << value << "\n";
		}

		std::cout << "\nEigenvectors:\n" << J;

		std::cout << "\nIterations: " << iter << "\n";

	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}


	return 0;
}