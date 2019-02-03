#ifndef GUARD_linear_element_h
#define GUARD_linear_element_h

#include <vector>

#include "element.h"

class LinearElement: public Element {
public:
    // Constructor
    using Element::Element;

    // Member functions
    double u(double) const;
    double previous_u(double) const;
    double d_u(double) const;
    double previous_d_u(double) const;

    // Basis functions
    double phi(double, unsigned) const;
    double d_phi(double, unsigned) const;
};

#endif
