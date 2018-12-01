#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), cos(), M_PI

#include "../src/Burgers_model.h"
#include "../src/utilities/utilities.h"

#include "projector_model.h"

double forcing_function(double x, double t)
{
    //return 1.0 + sin(1.0 * M_PI * x) * sin(1.0 * M_PI * t) +
    //             sin(2.0 * M_PI * x) * sin(2.0 * M_PI * t) +
    //             sin(3.0 * M_PI * x) * sin(3.0 * M_PI * t);

    //return sin(1.0 * M_PI * x) * sin(1.0 * M_PI * t) +
    //       sin(2.0 * M_PI * x) * sin(2.0 * M_PI * t) +
    //       sin(3.0 * M_PI * x) * sin(3.0 * M_PI * t);

    //return (1.0 + 2.0 * sin(1.5 * 2.0 * M_PI * x) * cos(1.5 * 2.0 * M_PI * t)) * cos(0.3 * 2.0 * M_PI * t);

    return 3.0 * sin(1.5 * 2.0 * M_PI * x) * cos(2.0 * M_PI * t);

    //return 1.0;
}

double initial_condition(double x)
{
    return 1.0 + sin(2.0 * M_PI * x);

    //return 1.0;
}

double left_boundary_value(double t)
{
    //return 1.0;

    return 1.0 + 0.75 * sin(0.2 * 2.0 * M_PI * t);
}

double right_boundary_value(double t)
{
    //return 1.0;

    return 1.0 - 0.75 * sin(0.2 * 2.0 * M_PI * t);
}

int main()
{
    // Problem parameters
    double nu            = 0.01;
    double x_left        = 0.0;
    double x_right       = 1.0;
    bool periodic_domain = false;
    double t_begin       = 0.0;
    double t_end         = 5.0;

    // Discretization
    unsigned N           = 1024;  // Number of elements
    unsigned timesteps   = 5001;
    unsigned N_proj      = 8;  // Number of element to project the solution onto
    unsigned proj_type   = 1;  // 1 for nodal projection, 2 for L2 projection

    // Define the x-coordinates of the nodes and the time steps
    std::vector<double> t = linspace(t_begin, t_end, timesteps);
    std::vector<double> x = linspace(x_left, x_right, N + 1);
    std::vector<double> x_proj = linspace(x_left, x_right, N_proj + 1);

    // Initialize the model
    BurgersModel model(N,
                       1,
                       nu,
                       x,
                       periodic_domain,
                       t_begin,
                       forcing_function,
                       initial_condition,
                       left_boundary_value,
                       right_boundary_value,
                       true);

    // Initialize the projection model
    ProjectorModel proj_model(model, N_proj, x_proj);

    // Allocate memory to store the solution
    std::vector<std::vector<double> > u(timesteps, std::vector<double>(N + 1));
    std::vector<std::vector<double> > u_proj(timesteps, std::vector<double>(N_proj + 1));

    // Compute and store the solution at t = 0
    u[0] = model.u(x);

    if (proj_type == 1) {
        u_proj[0] = model.u(x_proj);
    } else if (proj_type == 2) {
        proj_model.project();
        u_proj[0] = proj_model.u_proj(x_proj);
    }

    // Propagate through time and store the solutions at t > 0
    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);
        u[i] = model.u(x);

        // Project solution
        if (proj_type == 1) {
            // Nodal projection
            u_proj[i] = model.u(x_proj);
        } else if (proj_type == 2) {
            // L2 projection
            proj_model.project();
            u_proj[i] = proj_model.u_proj(x_proj);
        }
    }

    // Save results
    save_vector(t,      "dns_data/t_dns.dat");
    save_vector(x,      "dns_data/x_dns.dat");
    save_matrix(u,      "dns_data/u_dns.dat");
    save_vector(x_proj, "dns_data/x_proj.dat");
    save_matrix(u_proj, "dns_data/u_proj.dat");

    return 0;
}
