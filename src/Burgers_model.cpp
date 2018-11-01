#include <iostream> // std::cout, std::endl
#include <functional> // std::function
#include <vector> // std::vector
#include <array> // std::array
#include <cmath> // sqrt()
#include <numeric> // std::inner_product
#include <omp.h>
#include <algorithm>
#include <memory> // std::unique_ptr

#include "Burgers_model.h"
#include "elements/element.h"
#include "elements/linear_element.h"
#include "linear_algebra/tridiagonal_matrix.h"
#include "utilities/utilities.h"

using std::cout;
using std::endl;
using std::function;
using std::vector;
using std::array;
using std::inner_product;

// Constructor
BurgersModel::BurgersModel(unsigned number_of_elements,
                           unsigned polynomial_order,
                           double nu,
                           vector<double> nodes,
                           bool periodic_domain,
                           double time,
                           function<double(double, double)> forcing_function,
                           function<double(double)> initial_condition,
                           function<double(double)> left_boundary_value,
                           function<double(double)> right_boundary_value,
                           bool verbose)
    : number_of_elements(number_of_elements)
    , polynomial_order(polynomial_order)
    , nu(nu)
    , nodes(nodes)
    , periodic_domain(periodic_domain)
    , time(time)
    , forcing_function(forcing_function)
    , initial_condition(initial_condition)
    , left_boundary_value(left_boundary_value)
    , right_boundary_value(right_boundary_value)
    , verbose(verbose)
    , number_of_nodes(nodes.size())
    , nodes_per_element(polynomial_order + 1)
    , F(vector<double>(number_of_nodes))
    , J(TridiagonalMatrix(number_of_nodes))
{
    // Store `Element' objects in `elements'
    create_elements();

    project_initial_condition();

    // Set the number of OpenMP threads
    //omp_set_num_threads(1);
}

// double BurgersModel::forcing_function(double x, double t)
// {
//     return 0.0;
// }

// double BurgersModel::initial_condition(double x)
// {
//     return 1.0 + sin(2 * M_PI * x - 1.0);
// }

// double BurgersModel::left_boundary_value(double t)
// {
//     return 0.0;
// }

// double BurgersModel::right_boundary_value(double t)
// {
//     return 0.0;
// }

void BurgersModel::create_elements()
{
    for (unsigned i = 0; i != number_of_elements; ++i) {

        vector<double> local_nodes(nodes_per_element);
        vector<unsigned> indices(nodes_per_element);

        for (unsigned j = 0; j != nodes_per_element; ++j) {
            unsigned global_node_index = i * (nodes_per_element - 1) + j;

            local_nodes[j] = nodes[global_node_index];
            indices[j] = global_node_index;
        }

        elements.push_back(LinearElement(i, local_nodes, indices, coefficients, previous_coefficients));
    }
}

void BurgersModel::advance_in_time(double new_time)
{
    if (verbose) {
        cout << "Advancing from t = " << time << " to t = " << new_time << '.' << endl;
    }

    // Change timestamps
    previous_time = time;
    time = new_time;

    // Make a copy of `coefficients' for the computation of du_dt
    previous_coefficients = coefficients;

    if (!periodic_domain) {
        // Set boundary conditions
        coefficients.front() = left_boundary_value(time);
        coefficients.back() = right_boundary_value(time);
    }

    // Perform Newton iterations
    for (unsigned i = 0; i != 20; ++i) {
        // Assemble the residual vector `-F' and the Jacobian matrix `J'
        assemble_F();
        assemble_J();

        // Solve system
        if (periodic_domain) {
            J.solve_cyclic(F);
        } else {
            // Modify the system to enforce the boundary conditions
            constrain_system(J, F);

            J.solve(F); // Stores result in `F'
        }

        dxpy(F, coefficients); // Stores result in `coefficients'

        // Compute the 2-norm of `F'
        double norm = sqrt(inner_product(F.begin(), F.end(), F.begin(), 0.0));

        // Check for convergence
        if (norm < 1.0e-12) {
            if (verbose) {
                cout << "    Newton's method converged after " << i + 1 << " iterations." << endl;
            }

            // Terminate the loop
            break;
        }

        if (i == 19 && verbose) {
            cout << "    Newton's method FAILED to converge within " << i + 1 << " iterations." << endl;
        }
    }
}

void BurgersModel::assemble_F()
{
    std::fill(F.begin(), F.end(), 0.0);

    array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    array<double, 3> integration_points = {-1.0 * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};

    double dt = time - previous_time;

    // Loop over elements
    #pragma omp parallel for if (number_of_elements > 100)
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != nodes_per_element; ++i) {

            // Loop over integration points
            for (unsigned ip = 0; ip != 3; ++ip) {

                double integration_point = integration_points[ip];
                double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                // Evaluate at `x'
                double u     = element->u(x);
                double du_dt = (u - element->previous_u(x)) / dt;
                double du_dx = element->d_u(x);
                double phi   = element->phi(x, i);
                double d_phi = element->d_phi(x, i);
                double f     = forcing_function(x, time);

                double integrand = 0.0;

                // Term 1
                integrand += du_dt * phi;

                // Term 2
                integrand += -0.5 * u * u * d_phi;

                // Term 3
                integrand += nu * du_dx * d_phi;

                // Term 4
                integrand += -1.0 * f * phi;

                #pragma omp atomic
                F[element->indices[i]] -= weights[ip] * 0.5 * element->width * integrand;
            }
        }
    }
}

void BurgersModel::assemble_J()
{
    J.fill();

    array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    array<double, 3> integration_points = {-1.0 * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};

    double dt = time - previous_time;

    // Loop over elements
    #pragma omp parallel for if (number_of_elements > 100)
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != nodes_per_element; ++i) {

            // Loop over the local nodes in `element' a second time
            for (unsigned j = 0; j != nodes_per_element; ++j) {

                // Loop over integration points
                for (unsigned ip = 0; ip != 3; ++ip) {

                    double integration_point = integration_points[ip];
                    double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                    // Evaluate at `x'
                    double u       = element->u(x);
                    double phi_i   = element->phi(x, i);
                    double phi_j   = element->phi(x, j);
                    double d_phi_i = element->d_phi(x, i);
                    double d_phi_j = element->d_phi(x, j);

                    double integrand = 0.0;

                    // Term 1
                    integrand += 1.0 / dt * phi_j * phi_i;

                    // Term 2
                    integrand += -1.0 * u * phi_j * d_phi_i;

                    // Term 3
                    integrand += nu * d_phi_j * d_phi_i;

                    #pragma omp atomic
                    J.element(element->indices[i], element->indices[j]) += weights[ip] * 0.5 * element->width * integrand;
                }
            }
        }
    }
}

void BurgersModel::constrain_system(TridiagonalMatrix &mat,
                                    vector<double> &vec)
{
    mat.element(0, 0) = 1.0;
    mat.element(0, 1) = 0.0;
    mat.element(number_of_nodes - 1, number_of_nodes - 1) = 1.0;
    mat.element(number_of_nodes - 1, number_of_nodes - 2) = 0.0;

    vec.front() = 0.0;
    vec.back() = 0.0;
}

void BurgersModel::apply_boundary_conditions(TridiagonalMatrix &mat,
                                             vector<double>    &vec,
                                             double            value_left,
                                             double            value_right)
{
    mat.element(0, 0) = 1.0;
    mat.element(0, 1) = 0.0;
    mat.element(number_of_nodes - 1, number_of_nodes - 1) = 1.0;
    mat.element(number_of_nodes - 1, number_of_nodes - 2) = 0.0;

    vec.front() = value_left;
    vec.back() = value_right;
}

void BurgersModel::project_initial_condition()
{
    vector<double> b(number_of_nodes);
    TridiagonalMatrix M(number_of_nodes);

    // Assemble the right-hand side vector `b'
    assemble_b(b);

    // Assemble the mass matrix `M'
    assemble_M(M);

    if (periodic_domain) {
        // Solve system
        M.solve_cyclic(b);
    } else {
        // Apply boundary conditions
        apply_boundary_conditions(M, b, left_boundary_value(nodes.front()), right_boundary_value(nodes.back()));

        // Solve system
        M.solve(b);
    }

    // Set `coefficients' to `b'
    //coefficients = b;
    coefficients.swap(b);
}


void BurgersModel::assemble_b(vector<double>& b)
{
    std::fill(b.begin(), b.end(), 0.0);

    array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    array<double, 3> integration_points = {-1.0 * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};

    // Loop over elements
    #pragma omp parallel for if (number_of_elements > 100)
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != nodes_per_element; ++i) {

            // Loop over integration points
            for (unsigned ip = 0; ip != 3; ++ip) {

                double integration_point = integration_points[ip];
                double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                // Evaluate at `x'
                double phi_i = element->phi(x, i);
                double f     = initial_condition(x);

                double integrand = phi_i * f;

                #pragma omp atomic
                b[element->indices[i]] += weights[ip] * 0.5 * element->width * integrand;
            }
        }
    }
}

void BurgersModel::assemble_M(TridiagonalMatrix& M)
{
    array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    array<double, 3> integration_points = {-1.0 * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};

    M.fill();

    // Loop over elements
    #pragma omp parallel for if (number_of_elements > 100)
    for (auto element = elements.begin(); element < elements.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != nodes_per_element; ++i) {

            // Loop over the local nodes in `element' a second time
            for (unsigned j = 0; j != nodes_per_element; ++j) {

                // Loop over integration points
                for (unsigned ip = 0; ip != 3; ++ip) {

                    double integration_point = integration_points[ip];
                    double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                    // Evaluate at `x'
                    double phi_i = element->phi(x, i);
                    double phi_j = element->phi(x, j);

                    double integrand = phi_i * phi_j;

                    #pragma omp atomic
                    M.element(element->indices[i], element->indices[j]) += weights[ip] * 0.5 * element->width * integrand;
                }
            }
        }
    }
}

double BurgersModel::interpolate(double x)
{
    if (x < nodes.front() || x > nodes.back()) {
        return 0.0;
    } else if (x == nodes.back()) {
        return elements.back().u(x);
    } else {
        unsigned node_index = std::upper_bound(nodes.begin(), nodes.end(), x) - nodes.begin();
        unsigned element_index = node_index / polynomial_order - 1;

        return elements[element_index].u(x);
    }
}

vector<double> BurgersModel::interpolate(vector<double> &x)
{
    vector<double> result(x.size());

    for (unsigned i = 0; i < x.size(); ++i) {
        result[i] = interpolate(x[i]);
    }

    return result;
}

vector<double> BurgersModel::project(unsigned N)
{
    // Build mesh

    unsigned Nn = N + 1;

    std::vector<LinearElement> coarse_mesh;
    std::vector<double> x_nodes = linspace(nodes.front(), nodes.back(), Nn);

    std::vector<double> c(Nn);
    std::vector<double> prev_c(Nn);

    for (unsigned i = 0; i != N; ++i) {

        vector<double> local_x(2);
        vector<unsigned> indices(2);

        for (unsigned j = 0; j != 2; ++j) {
            local_x[j] = x_nodes[i + j];
            indices[j] = i + j;
        }

        coarse_mesh.push_back(LinearElement(i, local_x, indices, c, prev_c));
    }

    // Assemble `bb'

    vector<double> bb(Nn, 0.0);

    array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    array<double, 3> integration_points = {-1.0 * sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};

    // Loop over elements
    for (auto element = coarse_mesh.begin(); element < coarse_mesh.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != 2; ++i) {

            // Loop over integration points
            for (unsigned ip = 0; ip != 3; ++ip) {

                double integration_point = integration_points[ip];
                double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                // Evaluate at `x'
                double phi_i = element->phi(x, i);
                double f     = interpolate(x);

                double integrand = phi_i * f;

                bb[element->indices[i]] += weights[ip] * 0.5 * element->width * integrand;
            }
        }
    }

    // Assemble `MM'

    TridiagonalMatrix MM(Nn);
    MM.fill();

    // Loop over elements
    for (auto element = coarse_mesh.begin(); element < coarse_mesh.end(); ++element) {

        // Loop over the local nodes in `element'
        for (unsigned i = 0; i != 2; ++i) {

            // Loop over the local nodes in `element' a second time
            for (unsigned j = 0; j != 2; ++j) {

                // Loop over integration points
                for (unsigned ip = 0; ip != 3; ++ip) {

                    double integration_point = integration_points[ip];
                    double x = element->x_left + 0.5 * (1 + integration_point) * element->width;

                    // Evaluate at `x'
                    double phi_i = element->phi(x, i);
                    double phi_j = element->phi(x, j);

                    double integrand = phi_i * phi_j;

                    MM.element(element->indices[i], element->indices[j]) += weights[ip] * 0.5 * element->width * integrand;
                }
            }
        }
    }

    // Solve system

    if (periodic_domain) {
        // Solve system
        MM.solve_cyclic(bb);
    } else {
        // Apply boundary conditions
        apply_boundary_conditions(MM, bb, left_boundary_value(nodes.front()), right_boundary_value(nodes.back()));

        // Solve system
        MM.solve(bb);
    }

    return bb;
}
