#ifndef GUARD_Burgers_model_h
#define GUARD_Burgers_model_h

#include <functional> // std::function
#include <vector> // std::vector
#include <memory> // std::unique_ptr

#include "element.h"
#include "linear_element.h"
#include "linear_algebra.h"

class BurgersModel {
public:
    // Constructor
    BurgersModel(unsigned,
                 unsigned,
                 double,
                 std::vector<double>,
                 double,
                 bool=true);

    // Member functions

    static inline double forcing_function(double x, double t);

    static inline double initial_condition(double x);

    static inline double left_boundary_value(double t);

    static inline double right_boundary_value(double t);

    void create_elements(unsigned);

    void advance_in_time(double);

    void assemble_F();

    void assemble_J();

    double residual(double, double, double, auto&, unsigned);

    double jacobian(double, double, double, auto&, unsigned, unsigned);

    void constrain_system(TridiagonalMatrix&, std::vector<double>&);

    void apply_boundary_conditions(TridiagonalMatrix&, std::vector<double>&,
        double, double);

    void project_initial_condition();

    void assemble_b(std::vector<double>&);

    void assemble_M(TridiagonalMatrix&);

    double interpolate(double);

    std::vector<double> interpolate(std::vector<double>&);

    // Data members

    unsigned number_of_elements, polynomial_order;
    double nu;
    std::vector<double> nodes;
    double time;
    bool verbose;
    unsigned number_of_nodes, nodes_per_element;
    std::vector<double> F;
    TridiagonalMatrix J;

    std::vector<LinearElement> elements;
    std::vector<double> coefficients, previous_coefficients;
    double previous_time;
};

#endif
