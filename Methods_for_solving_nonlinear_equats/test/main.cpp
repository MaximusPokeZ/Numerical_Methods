#include "methods.h"

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


		long double eps;
		infile >> eps;

		infile.close();

		std::cout << "Dichotomy: First solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << dichotomy(eps, func8, 0.0, 1.0) << "\n";
		std::cout << "Dichotomy:Second solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << dichotomy(eps, func8, -0.99, 0.0) << "\n";

		std::cout << "Newton: First solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << newton_method(eps, 0.0, 1.0, func8, dfunc8, ddfunc8) << "\n";
		std::cout << "Newton: Second solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << newton_method(eps, -0.99, 0.0, func8, dfunc8, ddfunc8) << "\n";

		std::cout << "Secant: First solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << secant_method(eps, 0.0, 1.0, func8, ddfunc8) << "\n";
		std::cout << "Secant: Second solution of equation ln(x + 1) - 2 * x^2 + 1 is: " << secant_method(eps, -0.99, 0.0, func8, ddfunc8) << "\n";
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}


	return 0;
}