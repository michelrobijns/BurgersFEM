#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), M_PI
#include <iostream>

#include "Burgers_model.h"
#include "utilities/utilities.h"


double u_e(double x, double t)
{
    return sin(4.0 * M_PI * x) * sin(4.0 * M_PI * t);
}

std::vector<double> u_e(std::vector<double> x, double t)
{
    std::vector<double> result(x.size());

    for (unsigned i = 0; i != x.size(); ++i) {
        result[i] = u_e(x[i], t);
    }

    return result;
}

double u_t(double x, double t)
{
    return 4.0 * M_PI * sin(4.0 * M_PI * x) * cos(4.0 * M_PI * t);
}

double u_x(double x, double t)
{
    return 4.0 * M_PI * cos(4.0 * M_PI * x) * sin(4.0 * M_PI * t);
}

double u_xx(double x, double t)
{
    return -1.0 * 4.0 * M_PI * 4.0 * M_PI * sin(4.0 * M_PI * x) * sin(4.0 * M_PI * t);
}


double forcing_function(double x, double t)
{
    return u_t(x, t) + u_e(x, t) * u_x(x, t) - 0.01 * u_xx(x, t);
}

double initial_condition(double x)
{
    return 0.0;
}

double left_boundary_value(double t)
{
    return 0.0;
}

double right_boundary_value(double t)
{
    return 0.0;
}

int main()
{
    std::cout << "Hello, World!" << std::endl;

    // Declare parameters
    unsigned number_of_elements = 4096;
    unsigned polynomial_order   = 1;
    double nu                   = 0.01;
    double x_left               = 0.0;
    double x_right              = 1.0;
    bool periodic_domain        = false;
    double t_begin              = 0.0;
    double t_end                = 1.0;
    unsigned timesteps          = 10001;

    // Create a vector containing the x-coordinates of the nodes
    unsigned number_of_nodes = number_of_elements * polynomial_order + 1;
    std::vector<double> nodes = linspace(x_left, x_right, number_of_nodes);

    // Initialize the model
    BurgersModel model(number_of_elements,
                       polynomial_order,
                       nu,
                       nodes,
                       periodic_domain,
                       t_begin,
                       forcing_function,
                       initial_condition,
                       left_boundary_value,
                       right_boundary_value,
                       true);

    // Setup storage for the results
    std::vector<double> t = linspace(t_begin, t_end, timesteps);
    std::vector<double> x = linspace(x_left, x_right, 256);
    std::vector<std::vector<double> > u(timesteps, std::vector<double>(x.size()));
    std::vector<std::vector<double> > u_exact(timesteps, std::vector<double>(x.size()));

    // Store the present solution (the initial condition)
    u[0] = model.interpolate(x);
    u_exact[0] = u_e(x, t[0]);

    // Propagate through time and store the solutions
    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);

        u[i] = model.interpolate(x);
        u_exact[i] = u_e(x, t[i]);
    }

    save_vector(t, "test_data/t.dat");
    save_vector(x, "test_data/x.dat");
    save_matrix(u, "test_data/u.dat");
    save_matrix(u_exact, "test_data/u_e.dat");

    return 0;
}
