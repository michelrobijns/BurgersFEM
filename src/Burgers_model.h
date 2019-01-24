#ifndef GUARD_Burgers_model_h
#define GUARD_Burgers_model_h

#include <functional> // std::function
#include <vector>     // std::vector

#include "elements/element.h"
#include "elements/linear_element.h"
#include "linear_algebra/tridiagonal_matrix.h"

class BurgersModel {
public:
    // Constructor
    BurgersModel(unsigned                              number_of_elements,
                 unsigned                              polynomial_order,
                 double                                nu,
                 std::vector<double>                   nodes,
                 bool                                  periodic_domain,
                 double                                time,
                 std::function<double(double, double)> forcing_function,
                 std::function<double(double)>         initial_condition,
                 std::function<double(double)>         left_boundary_value,
                 std::function<double(double)>         right_boundary_value,
                 bool=true);

    // Member functions

    void create_elements();

    void advance_in_time(double new_time);

    void assemble_F();

    void assemble_J();

    void constrain_system(TridiagonalMatrix& mat, std::vector<double>& vec);

    void apply_boundary_conditions(TridiagonalMatrix& mat, std::vector<double>& vec,
        double value_left, double value_right);

    void project_initial_condition();

    void assemble_b(std::vector<double>& b);

    void assemble_M(TridiagonalMatrix& M);

    double u(double x);

    std::vector<double> u(std::vector<double>& x);

    double energy();

    // Data members

    unsigned number_of_elements, polynomial_order;
    double nu;
    std::vector<double> nodes;
    bool periodic_domain;
    double time;
    std::function<double(double, double)> forcing_function;
    std::function<double(double)> initial_condition;
    std::function<double(double)> left_boundary_value;
    std::function<double(double)> right_boundary_value;
    bool verbose;
    unsigned number_of_nodes, nodes_per_element;
    std::vector<double> F;
    TridiagonalMatrix J;

    std::vector<LinearElement> elements;
    std::vector<double> coefficients, previous_coefficients;
    double previous_time;
};

#endif
