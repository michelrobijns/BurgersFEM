#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), M_PI

#include "Burgers_model.h"
#include "utilities.h"

double forcing_function(double x, double t)
{
    return 0.0;
}

double initial_condition(double x)
{
    return 1.0 + sin(2 * M_PI * x - 1.0);
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
    // Declare parameters
    unsigned number_of_elements = 1024;
    unsigned polynomial_order   = 1;
    double nu                   = 0.03;
    double x_left               = 0.0;
    double x_right              = 1.0;
    bool periodic_domain        = true;
    double t_begin              = 0.0;
    double t_end                = 1.50;
    unsigned timesteps          = 151;

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
    std::vector<double> x = linspace(x_left, x_right, number_of_elements);
    std::vector<std::vector<double> > u(timesteps, std::vector<double>(number_of_elements));

    // Store the present solution (the initial condition)
    u[0] = model.interpolate(x);

    // Propagate through time and store the solutions
    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);

        u[i] = model.interpolate(x);
    }

    save_vector(t, "t.dat");
    save_vector(x, "x.dat");
    save_matrix(u, "u.dat");

    return 0;
}
