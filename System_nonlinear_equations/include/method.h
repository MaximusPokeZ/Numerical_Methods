#ifndef NUMERICAL_METHODS_METHOD_H
#define NUMERICAL_METHODS_METHOD_H

#include <cmath>
#include <functional>
#include "Matrix.h"

// first eq
long double f1 (long double x1, long double x2);
long double dx1_f1(long double x1, long double x2);
long double dx2_f1(long double x1, long double x2);

// second eq
long double f2 (long double x1, long double x2);
long double dx1_f2(long double x1, long double x2);
long double dx2_f2(long double x1, long double x2);


long double determinant_2x2(const Matrix<long double>& J);
size_t newton_method_system(const long double eps, long double& x1, long double& x2,
							Matrix<std::function<long double(long double, long double)>> J,
							Matrix<std::function<long double(long double, long double)>> A1,
							Matrix<std::function<long double(long double, long double)>> A2);


long double phi1 (long double x1, long double x2);
long double phi1_dx1 (long double x1, long double x2);
long double phi1_dx2 (long double x1, long double x2);
long double phi1_m (long double x1, long double x2);
long double phi1_dx1_m (long double x1, long double x2);
long double phi1_dx2_m (long double x1, long double x2);
long double phi2 (long double x1, long double x2);
long double phi2_dx1 (long double x1, long double x2);
long double phi2_dx2 (long double x1, long double x2);

size_t simple_iteration_method(long double eps, long double& x1, long double& x2,
							   std::function<long double(long double, long double)> phi1,
							   std::function<long double(long double, long double)> phi2, const long double& q);


#endif //NUMERICAL_METHODS_METHOD_H
