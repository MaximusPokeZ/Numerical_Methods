#include <iostream>
#include <fstream>
#include "gauss.h"

int main(int argc, char* argv[])
{

	if (argc != 2)
	{
		std::cerr << "No file provided";
		return -1;
	}

	std::ifstream in(argv[1]);

	if (!in.is_open())
	{
		std::cerr << "File not found";
		return -1;
	}

	size_t n;
	in >> n;

	Matrix<long double> A(n, n);
	std::vector<long double> b(n);

	in >> A;

	for (size_t i = 0; i < n; i++)
	{
		in >> b[i];
	}

	in.close();

	std::cout << "Matrix A:\n";
	std::cout << A;
	std::cout << "\nVector b:\n";
	for (const auto& val : b)
	{
		std::cout << val << "\n";
	}

	long double det = 1.0l;
	method_gauss(A, b, det);

	std::cout << "\nSolution:\n";
	for (const auto& val : b)
	{
		std::cout << val << "\n";
	}
	std::cout << "\nDeterminant: " << det << "\n";

	Matrix<long double> reverse = get_reverse_matrix(A);

	std::cout << "\nMatrix A^-1:\n";
	std::cout << reverse;


	std::cout << std::setprecision(-1) <<"\n";
	std::cout << A * reverse;

	return 0;
}
