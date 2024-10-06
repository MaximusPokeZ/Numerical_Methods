#include "method.h"

int main()
{
	try
	{
		//Newton

		Matrix<std::function<long double(long double, long double)>> J(2, 2);
		J[0][0] = dx1_f1;
		J[0][1] = dx2_f1;
		J[1][0] = dx1_f2;
		J[1][1] = dx2_f2;

		Matrix<std::function<long double(long double, long double)>> A1(2, 2);
		A1[0][0] = f1;
		A1[0][1] = dx2_f1;
		A1[1][0] = f2;
		A1[1][1] = dx2_f2;

		Matrix<std::function<long double(long double, long double)>> A2(2, 2);
		A2[0][0] = dx1_f1;
		A2[0][1] = f1;
		A2[1][0] = dx1_f2;
		A2[1][1] = f2;

		long double eps = 0.000001;


		long double a1 = -3, a2 = -1.5;
		long double b1 = -2.5, b2 = -1;
		long double x1 = (a1 + b1) * 0.5, x2 = (a2 + b2) * 0.5;

		size_t count_iter = newton_method_system(eps, x1, x2, J, A1, A2);

		std::cout << "Newton method first:\n";
		std::cout << "x1 = " << x1 << "\n";
		std::cout << "x2 = " << x2 << "\n";
		std::cout << "Number of iterations = " << count_iter << "\n\n";

		a1 = 2, a2 = 1.5;
		b1 = 2.5, b2 = 2;
		x1 = (a1 + b1) * 0.5, x2 = (a2 + b2) * 0.5;

		count_iter = newton_method_system(eps, x1, x2, J, A1, A2);

		std::cout << "Newton method second:\n";
		std::cout << "x1 = " << x1 << "\n";
		std::cout << "x2 = " << x2 << "\n";
		std::cout << "Number of iterations = " << count_iter << "\n\n";

		// Simple iter
		a1 = 2, a2 = 1.5;
		b1 = 2.5, b2 = 1.7;
		x1 = (a1 + b1) * 0.5, x2 = (a2 + b2) * 0.5;

		Matrix<std::function<long double(long double, long double)>> phi_mat(2, 2);
		phi_mat[0][0] = phi1_dx1;
		phi_mat[0][1] = phi1_dx2;
		phi_mat[1][0] = phi2_dx1;
		phi_mat[1][1] = phi2_dx2;

		Matrix<long double> phi_val(2, 2);
		phi_val[0][0] = phi_mat[0][0](b1, b2); phi_val[0][1] = phi_mat[0][1](b1, b2);
		phi_val[1][0] = phi_mat[1][0](b1, b2); phi_val[1][1] = phi_mat[1][1](b1, b2);

		long double q = std::abs(phi_val.row_norm());
		if (q >= 1) throw std::runtime_error("The convergence condition is not satisfied!");

		count_iter = simple_iteration_method(eps, x1, x2, phi1, phi2, q);

		std::cout << "Simple iter method first:\n";
		std::cout << "x1 = " << x1 << "\n";
		std::cout << "x2 = " << x2 << "\n";
		std::cout << "Number of iterations = " << count_iter << "\n\n";


		a1 = -3, a2 = -1.5;
		b1 = -2.5, b2 = -1;
		x1 = (a1 + b1) * 0.5, x2 = (a2 + b2) * 0.5;

		phi_mat[0][0] = phi1_dx1_m;
		phi_mat[0][1] = phi1_dx2_m;
		phi_val[0][0] = phi_mat[0][0](b1, b2); phi_val[0][1] = phi_mat[0][1](b1, b2);
		phi_val[1][0] = phi_mat[1][0](b1, b2); phi_val[1][1] = phi_mat[1][1](b1, b2);
		q = std::abs(phi_val.row_norm());
		if (q >= 1) throw std::runtime_error("The convergence condition is not satisfied!");

		count_iter = simple_iteration_method(eps, x1, x2, phi1_m, phi2, q);


		std::cout << "Simple iter method second:\n";
		std::cout << "x1 = " << x1 << "\n";
		std::cout << "x2 = " << x2 << "\n";
		std::cout << "Number of iterations = " << count_iter;








	} catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}










	return 0;
}