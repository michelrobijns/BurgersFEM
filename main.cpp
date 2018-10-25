#include <vector> // std::vector
#include <cmath> // sqrt(), sin(), cos()

#include "linear_element.h"
#include "Burgers_model.h"
#include "utilities.h"

using std::vector;


int main()
{
    // Declare parameters
    unsigned number_of_elements = 1024;
    unsigned polynomial_order   = 1;
    double nu                   = 0.03;
    double x_left               = 0.0;
    double x_right              = 1.0;
    double t_begin              = 0.0;
    double t_end                = 1.50;
    double timestep             = 0.01;

    // Create a vector containing the x-coordinates of the nodes
    unsigned number_of_nodes = number_of_elements * polynomial_order + 1;
    vector<double> nodes = linspace(x_left, x_right, number_of_nodes);

    // Initialize the model
    BurgersModel model(number_of_elements, polynomial_order, nu, nodes,
        t_begin, false);

    // Setup storage for the results
    vector<double> time = {t_begin};
    vector<double> x = linspace(x_left, x_right, 256);

    // Store the present solution (the initial condition)
    vector<vector<double> > u = {model.interpolate(x)};

    // Advance the model in time and the store the solutions
    for (double t = t_begin + timestep; t < t_end + timestep; t += timestep) {
        model.advance_in_time(t);

        time.push_back(t);
        u.push_back(model.interpolate(x));  // Note to self: remove `push_back' because the final size is known at compile time
    }

    save_vector(time, "time.dat");
    save_vector(x, "x.dat");
    save_matrix(u, "u.dat");

    return 0;
}
