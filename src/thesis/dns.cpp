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
    unsigned N_write        = 1024;

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

    // Storage
    std::vector<double> t_write = linspace(t_begin, t_end,   timesteps);
    std::vector<double> x_write = linspace(x_left,  x_right, N_write + 1);
    std::vector<double> x_proj  = linspace(x_left,  x_right, N_proj + 1);

    std::vector<std::vector<double> > u        (timesteps, std::vector<double>(N_write + 1));
    std::vector<std::vector<double> > u_proj   (timesteps, std::vector<double>(N_proj  + 1));
    std::vector<std::vector<double> > u_prime  (timesteps, std::vector<double>(N_write + 1));
    std::vector<std::vector<double> > u_prime_t(timesteps, std::vector<double>(N_write + 1));

    // Initialize the projection model
    ProjectorModel proj_model(model, N_proj, x_proj);


    // Compute and store the solution at t = 0

    proj_model.project();

    for (unsigned j = 0; j < (N_write + 1); j++) {
        double xj = x_write[j];

        double u_xj      = model.u(xj);
        double u_proj_xj = proj_model.u_proj(xj);

        u[0][j]         = u_xj;
        u_prime[0][j]   = u_xj - u_proj_xj;
        u_prime_t[0][j] = 0.0;
    }

    u_proj[0] = proj_model.u_proj(x_proj);

    // Propagate through time and store the solutions at t > 0

    for (unsigned i = 1; i != timesteps; i++) {
        model.advance_in_time(t[i]);

        // Store solution

        proj_model.project();

        // Compute u' and u'_t

        for (unsigned j = 0; j < (N_write + 1); j++) {
            double xj = x_write[j];

            double u_xj       = model.u(xj);
            double u_proj_xj  = proj_model.u_proj(xj);

            u[i][j]         = u_xj;
            u_prime[i][j]   = u_xj - u_proj_xj;
            u_prime_t[i][j] = (u_prime[i][j] - u_prime[i - 1][j]) /
                              (t[i] - t[i - 1]);
        }

        u_proj[i] = proj_model.u_proj(x_proj);
    }

    // Save results

    save_vector(t_write, "dns_data/t.dat");
    save_vector(x_write, "dns_data/x.dat");
    save_matrix(u, "dns_data/u.dat");
    save_vector(x_proj, "dns_data/x_proj.dat");
    save_matrix(u_proj, "dns_data/u_proj.dat");
    save_matrix(u_prime, "dns_data/u_prime.dat");
    save_matrix(u_prime_t, "dns_data/u_prime_t.dat");

    return 0;
}
