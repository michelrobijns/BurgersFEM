#include <vector>  // std::vector
#include <algorithm>  // std::find

#include "basis_functions.h"

using std::vector;
using std::find;

double piecewise_linear(double x, unsigned i, vector<double> nodes)
{
    if (i == 0) {
        return (x - nodes[1]) / (nodes[0] - nodes[1]);
    } else {
        return (x - nodes[0]) / (nodes[1] - nodes[0]);
    }
}

vector<double> piecewise_linear(vector<double> x, unsigned i, vector<double> nodes)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = piecewise_linear(x[j], i, nodes);
    }

    return results;
}

double d_piecewise_linear(double x, unsigned i, vector<double> nodes)
{
    if (i == 0) {
        return 1.0 / (nodes[0] - nodes[1]);
    } else {
        return 1.0 / (nodes[1] - nodes[0]);
    }
}

vector<double> d_piecewise_linear(vector<double> x, unsigned i, vector<double> nodes)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = d_piecewise_linear(x[j], i, nodes);
    }

    return results;
}

double d_piecewise_linear(double x, unsigned i, vector<double> nodes, unsigned order)
{
    if (order > 1) {
        return 0.0;
    } else if (i == 0) {
        return 1.0 / (nodes[0] - nodes[1]);
    } else {
        return 1.0 / (nodes[1] - nodes[0]);
    }
}

vector<double> d_piecewise_linear(vector<double> x, unsigned i, vector<double> nodes, unsigned order)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = d_piecewise_linear(x[j], i, nodes, order);
    }

    return results;
}

double lagrange(double x, unsigned i, vector<double> nodes)
{
    typedef vector<double>::size_type vec_sz;
    vec_sz size = nodes.size();

    double result = 1.0;

    for (unsigned j = 0; j != size; ++j) {
        if (j != i) {
            result *= ((x - nodes[j]) / (nodes[i] - nodes[j]));
        }
    }

    return result;
}

vector<double> lagrange(vector<double> x, unsigned i, vector<double> nodes)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = lagrange(x[j], i, nodes);
    }

    return results;
}

double d_lagrange(double x, unsigned i, vector<double> nodes)
{
    typedef vector<double>::size_type vec_sz;
    vec_sz size = nodes.size();

    double result = 0.0;

    for (unsigned j = 0; j != size; ++j) {
        if (j != i) {
            double partial_result = 1.0 / (nodes[i] - nodes[j]);

            for (unsigned k = 0; k != size; ++k) {
                if (k != i && k != j) {
                    partial_result *= ((x - nodes[k]) / (nodes[i] - nodes[k]));
                }
            }

            result += partial_result;
        }
    }

    return result;
}

vector<double> d_lagrange(vector<double> x, unsigned i, vector<double> nodes)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = d_lagrange(x[j], i, nodes);
    }

    return results;
}

double d_lagrange(double x, unsigned i, vector<double> nodes, unsigned order,
    vector<unsigned> nodes_to_exclude)
{
    typedef vector<double>::size_type vec_sz;
    vec_sz size = nodes.size();

    double result = 0.0;

    for (unsigned j = 0; j != size; ++j) {
        if (j != i && !contains(nodes_to_exclude, j)) {
            double partial_result = 1.0 / (nodes[i] - nodes[j]);

            if (order != 1) {
                // Make a copy of `nodes_to_exclude' before adding the new node
                // that must be excluded.  Copying `nodes_to_exclude' first is
                // necessary because the same vector (as in, the same memory) is
                // used in subsequent iterations of the for loop.
                vector<unsigned> nodes_to_exclude_copy(nodes_to_exclude);
                nodes_to_exclude_copy.push_back(j);

                result += partial_result * d_lagrange(x, i, nodes, order - 1,
                    nodes_to_exclude_copy);
            } else {
                for (unsigned k = 0; k != size; ++k) {
                    if (k != i && k != j && !contains(nodes_to_exclude, k)) {
                        partial_result *= ((x - nodes[k]) / (nodes[i] - nodes[k]));
                    }
                }

                result += partial_result;
            }
        }
    }

    return result;
}

vector<double> d_lagrange(vector<double> x, unsigned i, vector<double> nodes, unsigned order,
    vector<unsigned> nodes_to_exclude)
{
    vector<double> results(x.size());

    for (vector<double>::size_type j = 0; j != x.size(); ++j) {
        results[j] = d_lagrange(x[j], i, nodes, order);
    }

    return results;
}

inline bool contains(vector<unsigned> vec, unsigned element)
{
    return find(vec.begin(), vec.end(), element) != vec.end();
}
