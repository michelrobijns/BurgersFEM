#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), M_PI

#include "../Burgers_model.h"
#include "../utilities/utilities.h"

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
    unsigned number_of_elements = 4096;
    unsigned polynomial_order   = 1;
    double nu                   = 0.01;
    double x_left               = 0.0;
    double x_right              = 1.0;
    bool periodic_domain        = true;
    double t_begin              = 0.0;
    double t_end                = 1.50;
    unsigned timesteps          = 1501;

    // Only store a subset of the solution
    unsigned storage_nodes      = 256;
    unsigned skip_timesteps     = 10;

    // Create a vector containing the x-coordinates of the nodes and the time
    // steps
    unsigned number_of_nodes = number_of_elements * polynomial_order + 1;
    std::vector<double> nodes = linspace(x_left, x_right, number_of_nodes);
    std::vector<double> time = linspace(t_begin, t_end, timesteps);

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
    unsigned storage_timesteps = (timesteps - 1) / skip_timesteps + 1;

    std::vector<double> t = linspace(t_begin, t_end, storage_timesteps);
    std::vector<double> x = linspace(x_left, x_right, storage_nodes);
    std::vector<std::vector<double> > u(storage_timesteps, std::vector<double>(storage_nodes));

    std::vector<double> x_proj = linspace(x_left, x_right, 9);
    std::vector<std::vector<double> > u_proj(storage_timesteps, std::vector<double>(9));

    // Store the present solution (the initial condition)
    u[0] = model.interpolate(x);
    u_proj[0] = model.project(8);

    // Propagate through time and store the solutions
    unsigned j = 1;

    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(time[i]);

        if (i % skip_timesteps == 0) {
            u[j] = model.interpolate(x);
            u_proj[j++] = model.project(8);
        }
    }

    // Save results to file
    save_vector(t, "dns_data/t.dat");
    save_vector(x, "dns_data/x.dat");
    save_matrix(u, "dns_data/u.dat");
    save_vector(x_proj, "dns_data/x_proj.dat");
    save_matrix(u_proj, "dns_data/u_proj.dat");

    return 0;
}
