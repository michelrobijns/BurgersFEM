#include <vector> // std::vector
#include <algorithm> // std::fill, std::upper_bound

#include "projector_model.h"
#include "../utilities/legendre_rule.h"

// Constructor
ProjectorModel::ProjectorModel(BurgersModel& model,
                               unsigned number_of_elements,
                               std::vector<double> nodes)
    : model(model)
    , number_of_elements(number_of_elements)
    , nodes(nodes)
    , number_of_nodes(number_of_elements + 1)
    , b(std::vector<double>(number_of_nodes))
    , A(TridiagonalMatrix(number_of_nodes))
{
    create_elements();
}

void ProjectorModel::create_elements()
{
    for (unsigned i = 0; i != number_of_elements; ++i) {

        std::vector<double> local_nodes(2);
        std::vector<unsigned> indices(2);

        for (unsigned j = 0; j != 2; ++j) {
            local_nodes[j] = nodes[i + j];
            indices[j] = i + j;
        }

        elements.push_back(LinearElement(i, local_nodes, indices, coefficients, previous_coefficients));
    }
}

void ProjectorModel::project()
{
    assemble_b();

    assemble_A();

    if (model.periodic_domain) {
        A.solve_cyclic(b);
    } else {
        apply_boundary_conditions(A, b, model.u(nodes.front()),
            model.u(nodes.back()));

        A.solve(b);
    }

    // We can efficiently swap `coefficients' with `b' because we will
    // reassemble `b' the next time `project' is called
    //coefficients.swap(b);
    coefficients = b;
}

void ProjectorModel::assemble_b()
{
    unsigned num_int_pts = 3;
    std::vector<double> int_weights;
    std::vector<double> int_pts;
    legendre_rule(num_int_pts, int_weights, int_pts);

    std::fill(b.begin(), b.end(), 0.0);

    // Loop over elements
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != 2; ++i) {

            // Loop over integration points
            for (unsigned k = 0; k != num_int_pts; ++k) {

                double x_ip = int_pts[k];
                double x = element->x_left + 0.5 * (1 + x_ip) * element->width;

                // Evaluate at `x'
                double phi_i = element->phi(x, i);
                double f     = model.u(x);

                double integrand = phi_i * f;

                b[element->indices[i]] +=
                    int_weights[k] * 0.5 * element->width * integrand;
            }
        }
    }
}

void ProjectorModel::assemble_A()
{
    unsigned num_int_pts = 3;
    std::vector<double> int_weights;
    std::vector<double> int_pts;
    legendre_rule(num_int_pts, int_weights, int_pts);

    A.fill();

    // Loop over elements
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != 2; ++i) {

            // Loop over the local nodes in `element' a second time
            for (unsigned j = 0; j != 2; ++j) {

                // Loop over integration points
                for (unsigned k = 0; k != num_int_pts; ++k) {

                    double x_ip = int_pts[k];
                    double x = element->x_left +
                        0.5 * (1 + x_ip) * element->width;

                    // Evaluate at `x'
                    double phi_i = element->phi(x, i);
                    double phi_j = element->phi(x, j);

                    double integrand = phi_i * phi_j;

                    A.element(element->indices[i], element->indices[j]) +=
                        int_weights[k] * 0.5 * element->width * integrand;
                }
            }
        }
    }
}

double ProjectorModel::u_proj(double x)
{
    if (x < nodes.front() || x > nodes.back()) {
        return 0.0;
    } else if (x == nodes.back()) {
        return elements.back().u(x);
    } else {
        unsigned node_index = std::upper_bound(nodes.begin(), nodes.end(), x) - nodes.begin();
        unsigned element_index = node_index - 1;

        return elements[element_index].u(x);
    }
}

std::vector<double> ProjectorModel::u_proj(std::vector<double> &x)
{
    std::vector<double> result(x.size());

    for (unsigned i = 0; i < x.size(); ++i) {
        result[i] = u_proj(x[i]);
    }

    return result;
}

void ProjectorModel::apply_boundary_conditions(TridiagonalMatrix &matrix,
    std::vector<double> &vector, double boundary_value_left,
    double boundary_value_right)
{
    matrix.element(0, 0) = 1.0;
    matrix.element(0, 1) = 0.0;

    matrix.element(number_of_nodes - 1, number_of_nodes - 1) = 1.0;
    matrix.element(number_of_nodes - 1, number_of_nodes - 2) = 0.0;

    vector.front() = boundary_value_left;
    vector.back() = boundary_value_right;
}
