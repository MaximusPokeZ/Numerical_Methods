//
// Created by Vladimir Zaslavtsev on 04.10.2024.
//

#ifndef NUMERICAL_METHODS_METHODS_H
#define NUMERICAL_METHODS_METHODS_H

#include <cmath>
#include <functional>

typedef long double (*func)(long double);

long double dichotomy (long double eps, func f, long double a, long double b);

long double newton_method(long double eps, long double a, long double b, std::function<long double(long double)> func,
						  std::function<long double(long double)> dfunc, std::function<long double(long double)> ddfunc);

long double secant_method(long double eps, long double a, long double b, std::function<long double(long double)> func, std::function<long double(long double)> ddfunc);

long double simple_iter (const std::function<long double(long double)>& phi_f, const std::function<long double(long double)>& d_phi_f, long double a, long double b, long double eps);

long double func8 (long double x);
long double dfunc8(long double x);
long double ddfunc8(long double x);
long double phi_x_p(long double x);
long double dphi_x_p (long double x);

#endif //NUMERICAL_METHODS_METHODS_H
