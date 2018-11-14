#ifndef GUARD_projector_model_h
#define GUARD_projector_model_h

#include <vector> // std::vector

#include "../src/elements/linear_element.h"
#include "../src/linear_algebra/tridiagonal_matrix.h"
#include "../src/Burgers_model.h"

class ProjectorModel {
public:
    // Constructor
    ProjectorModel(BurgersModel& model,
                   unsigned number_of_elements,
                   std::vector<double> nodes);

    // Member functions

    void create_elements();

    void project();

    void assemble_b();

    void assemble_A();

    double u_proj(double x);

    std::vector<double> u_proj(std::vector<double>& x);

    void apply_boundary_conditions(TridiagonalMatrix& matrix,
        std::vector<double>& vector, double left_boundary_value,
        double right_boundary_value);

    // Data members

    BurgersModel& model;
    unsigned number_of_elements;
    std::vector<double> nodes;

    unsigned number_of_nodes;
    std::vector<double> b;
    TridiagonalMatrix A;
    std::vector<LinearElement> elements;
    std::vector<double> coefficients;
    std::vector<double> previous_coefficients;
};

#endif
