#ifndef GUARD_utilities_h
#define GUARD_utilities_h

#include <vector>
#include <string>

void print_vector(std::vector<double>, unsigned=3, unsigned=3, unsigned=100, bool=true);

template <class T>
void print_matrix(T, unsigned=3, unsigned=3, unsigned=100, bool=true);

void save_vector(std::vector<double>, std::string);

void save_matrix(std::vector<std::vector<double> >, std::string, std::string=" ");

std::vector<double> linspace(double, double, unsigned, bool=true);

#endif
