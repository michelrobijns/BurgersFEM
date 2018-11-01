#include <vector>

#include "element.h"
#include "basis_functions.h"

using std::vector;

// Constructor
Element::Element(unsigned index,
                 vector<double> nodes,
                 vector<unsigned> indices,
                 vector<double> &coefficients,
                 vector<double> &previous_coefficients)
    : index(index)
    , nodes(nodes)
    , indices(indices)
    , coefficients(coefficients)
    , previous_coefficients(previous_coefficients)
{
    x_left = nodes.front();
    x_right = nodes.back();

    width = x_right - x_left;
}

// Member functions

double Element::u(double x) const // Why does this not work without the `const' keyword?
{
    double result = 0.0;

    // Loop over the local indices of the local nodes
    for (unsigned i = 0; i != nodes.size(); ++i) {
        result += coefficients[indices[i]] * phi(x, i);
    }

    return result;
}

double Element::previous_u(double x) const // Why does this not work without the `const' keyword?
{
    double result = 0.0;

    // Loop over the local indices of the local nodes
    for (unsigned i = 0; i != nodes.size(); ++i) {
        result += previous_coefficients[indices[i]] * phi(x, i);
    }

    return result;
}

double Element::d_u(double x) const
{
    double result = 0.0;

    // Loop over the local indices of the local nodes
    for (unsigned i = 0; i != nodes.size(); ++i) {
        result += coefficients[indices[i]] * d_phi(x, i);
    }

    return result;
}

double Element::d_u(double x, unsigned order) const
{
    double result = 0.0;

    // Loop over the local indices of the local nodes
    for (unsigned i = 0; i != nodes.size(); ++i) {
        result += coefficients[indices[i]] * d_phi(x, i, order);
    }

    return result;
}

// Basis functions

double Element::phi(double x, unsigned i) const
{
    //return lagrange(x, i, nodes);
    return piecewise_linear(x, i, nodes);
}

vector<double> Element::phi(vector<double> x, unsigned i) const
{
    //return lagrange(x, i, nodes);
    return piecewise_linear(x, i, nodes);
}

double Element::d_phi(double x, unsigned i) const
{
    //return d_lagrange(x, i, nodes);
    return d_piecewise_linear(x, i, nodes);
}

vector<double> Element::d_phi(vector<double> x, unsigned i) const
{
    //return d_lagrange(x, i, nodes);
    return d_piecewise_linear(x, i, nodes);
}

double Element::d_phi(double x, unsigned i, unsigned order) const
{
    //return d_lagrange(x, i, nodes, order);
    return d_piecewise_linear(x, i, nodes, order);
}

vector<double> Element::d_phi(vector<double> x, unsigned i, unsigned order) const
{
    //return d_lagrange(x, i, nodes, order);
    return d_piecewise_linear(x, i, nodes, order);
}
