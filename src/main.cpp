#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), M_PI

#include "Burgers_model.h"
#include "utilities/utilities.h"

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
    double nu                   = 0.01;
    double x_left               = 0.0;
    double x_right              = 1.0;
    bool periodic_domain        = true;
    double t_begin              = 0.0;
    double t_end                = 1.50;
    unsigned timesteps          = 151;

    // Only store a subset of the solution
    unsigned storage_nodes      = 256;
    unsigned skip_timesteps     = 1;

    // Create a vector containing the x-coordinates of the nodes and the time
    // steps
    unsigned number_of_nodes = number_of_elements * polynomial_order + 1;
    std::vector<double> x = linspace(x_left, x_right, number_of_nodes);
    std::vector<double> t = linspace(t_begin, t_end, timesteps);

    // Initialize the model
    BurgersModel model(number_of_elements,
                       polynomial_order,
                       nu,
                       x,
                       periodic_domain,
                       t_begin,
                       forcing_function,
                       initial_condition,
                       left_boundary_value,
                       right_boundary_value,
                       false);

    // Setup storage for the results
    unsigned storage_timesteps = (timesteps - 1) / skip_timesteps + 1;

    std::vector<double> t_storage = linspace(t_begin, t_end, storage_timesteps);
    std::vector<double> x_storage = linspace(x_left, x_right, storage_nodes);
    std::vector<std::vector<double> > u_storage(storage_timesteps, std::vector<double>(storage_nodes));

    // Store the present solution (the initial condition)
    u_storage[0] = model.u(x_storage);

    // Propagate through time and store the solutions
    unsigned j = 1;

    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);

        if (i % skip_timesteps == 0) {
            u_storage[j++] = model.u(x_storage);
        }
    }

    // Save results to file
    save_vector(t_storage, "data/t.dat");
    save_vector(x_storage, "data/x.dat");
    save_matrix(u_storage, "data/u.dat");

    return 0;
}
