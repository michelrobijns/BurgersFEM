#include <vector>
#include <iostream>

#include "linear_element.h"
#include "basis_functions.h"

using std::vector;

// Member functions

double LinearElement::u(double x) const // Why does this not work without the `const' keyword?
{
    return coefficients[indices[0]] * phi(x, 0) + coefficients[indices[1]] * phi(x, 1);
}

double LinearElement::previous_u(double x) const // Why does this not work without the `const' keyword?
{
    return previous_coefficients[indices[0]] * phi(x, 0) + previous_coefficients[indices[1]] * phi(x, 1);
}

double LinearElement::d_u(double x) const
{
    return coefficients[indices[0]] * d_phi(x, 0) + coefficients[indices[1]] * d_phi(x, 1);
}

// Basis functions

double LinearElement::phi(double x, unsigned i) const
{
    // This function was manually inlined from basis_functions.h
    if (i == 0) {
        return (x - nodes[1]) / (nodes[0] - nodes[1]);
    } else {
        return (x - nodes[0]) / (nodes[1] - nodes[0]);
    }
}

double LinearElement::d_phi(double x, unsigned i) const
{
    // This function was manually inlined from basis_functions.h
    if (i == 0) {
        return 1.0 / (nodes[0] - nodes[1]);
    } else {
        return 1.0 / (nodes[1] - nodes[0]);
    }
}
