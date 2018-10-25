#ifndef GUARD_element_h
#define GUARD_element_h

#include <vector>

class Element {
public:
    // Constructor
    Element(unsigned,
            std::vector<double>,
            std::vector<unsigned>,
            std::vector<double>&,
            std::vector<double>&);

    // Member functions
    double u(double) const;
    double previous_u(double) const;
    double d_u(double) const;
    double d_u(double, unsigned) const;

    // Basis functions
    double phi(double, unsigned) const;
    std::vector<double> phi(std::vector<double>, unsigned) const;
    double d_phi(double, unsigned) const;
    std::vector<double> d_phi(std::vector<double>, unsigned) const;
    double d_phi(double, unsigned, unsigned) const;
    std::vector<double> d_phi(std::vector<double>, unsigned, unsigned) const;

    // Data members
    unsigned index;
    double x_left, x_right, width;
    std::vector<double> nodes;
    std::vector<unsigned> indices;
    std::vector<double> &coefficients, &previous_coefficients;
};

#endif
