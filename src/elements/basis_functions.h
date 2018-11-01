#ifndef GUARD_basis_functions_h
#define GUARD_basis_functions_h

#include <vector>

double piecewise_linear(double, unsigned, std::vector<double>);

std::vector<double> piecewise_linear(std::vector<double>, unsigned, std::vector<double>);

double d_piecewise_linear(double, unsigned, std::vector<double>);

std::vector<double> d_piecewise_linear(std::vector<double>, unsigned, std::vector<double>);

double d_piecewise_linear(double, unsigned, std::vector<double>, unsigned);

std::vector<double> d_piecewise_linear(std::vector<double>, unsigned, std::vector<double>, unsigned);

double lagrange(double, unsigned, std::vector<double>);

std::vector<double> lagrange(std::vector<double>, unsigned, std::vector<double>);

double d_lagrange(double, unsigned, std::vector<double>);

std::vector<double> d_lagrange(std::vector<double>, unsigned, std::vector<double>);

double d_lagrange(double, unsigned, std::vector<double>, unsigned, std::vector<unsigned> = std::vector<unsigned>());

std::vector<double> d_lagrange(std::vector<double>, unsigned, std::vector<double>, unsigned, std::vector<unsigned> = std::vector<unsigned>());

bool contains(std::vector<unsigned>, unsigned);

#endif
