#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), M_PI

#include "../Burgers_model.h"
#include "../utilities/utilities.h"

#include "projector_model.h"

double forcing_function(double x, double t)
{
    return 0.0;
}

double initial_condition(double x)
{
    return 1.0 + sin(2.0 * M_PI * x - 1.0);
}

double left_boundary_value(double t)
{
    return 1.0;
}

double right_boundary_value(double t)
{
    return 1.0;
}

int main()
{
    // Problem parameters
    double nu               = 0.01;
    double x_left           = 0.0;
    double x_right          = 1.0;
    bool periodic_domain    = true;
    double t_begin          = 0.0;
    double t_end            = 2.0;

    // Discretization
    unsigned N              = 1024;  // Number of elements
    unsigned timesteps      = 2001;
    unsigned N_proj         = 8;  // Number of element to project the solution onto

    // Parameters of the solution that is written to disk
    unsigned N_write        = 256;
    unsigned write_interval = 10;  // Number of timesteps between solution writes

    // Define the x-coordinates of the nodes and the time steps
    std::vector<double> x = linspace(x_left, x_right, N + 1);
    std::vector<double> t = linspace(t_begin, t_end, timesteps);

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
                       false);

    // Solution
    unsigned write_timesteps = (timesteps - 1) / write_interval + 1;

    std::vector<double> t_write = linspace(t_begin, t_end, write_timesteps);
    std::vector<double> x_write = linspace(x_left, x_right, N_write + 1);

    std::vector<std::vector<double>> u_write(write_timesteps,
        std::vector<double>(N_write + 1));

    // Projected (large-scale) solution
    std::vector<double> x_proj = linspace(x_left, x_right, N_proj + 1);
    std::vector<std::vector<double> > u_proj(write_timesteps,
        std::vector<double>(N_proj + 1));

    // Small-scale solution
    std::vector<std::vector<double> > u_prime(write_timesteps,
        std::vector<double>(N_write + 1));
    std::vector<std::vector<double> > u_prime_t(write_timesteps,
        std::vector<double>(N_write + 1));

    // Initialize the projection model
    ProjectorModel proj_model(model, N_proj, x_proj);






    // Store the present solution (the initial condition)
    u_write[0] = model.u(x_write);

    proj_model.project();

    u_proj[0] = proj_model.u_proj(x_proj);

    // Compute u_prime and u_prime_t
    for (unsigned k = 0; k < (N_write + 1); k++) {
        double xk = x_write[k];

        u_prime[0][k]   = model.u(xk) - proj_model.u_proj(xk);
        u_prime_t[0][k] = 0.0;
    }






    // Propagate through time and store the solutions
    unsigned j = 1;

    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);

        if ((i + 1) % write_interval == 0) {
            proj_model.project();

            // Compute u_prime
            for (unsigned k = 0; k < (N_write + 1); k++) {
                double xk = x_write[k];

                u_prime[j][k] = model.u(xk) - proj_model.u_proj(xk);
            }
        }



        if (i % write_interval == 0) {
            u_write[j] = model.u(x_write);

            proj_model.project();

            u_proj[j] = proj_model.u_proj(x_proj);

            // Compute uprime and uprime_t
            for (unsigned k = 0; k < (N_write + 1); k++) {
                double xk = x_write[k];

                double u_prime_temp = model.u(xk) - proj_model.u_proj(xk);

                u_prime_t[j][k] = (u_prime_temp - u_prime[j][k]) /
                                  (t[i] - t[i-1]);

                u_prime[j][k] = u_prime_temp;
            }

            // Increment `j'
            j++;
        }
    }

    // Save results
    save_vector(t_write, "dns_data/t.dat");
    save_vector(x_write, "dns_data/x.dat");
    save_matrix(u_write, "dns_data/u.dat");
    save_vector(x_proj, "dns_data/x_proj.dat");
    save_matrix(u_proj, "dns_data/u_proj.dat");
    save_matrix(u_prime, "dns_data/u_prime.dat");
    save_matrix(u_prime_t, "dns_data/u_prime_t.dat");

    return 0;
}
