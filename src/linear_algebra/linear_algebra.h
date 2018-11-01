#ifndef GUARD_linear_algebra_h
#define GUARD_linear_algebra_h

#include <vector>
#include <array>

// class TridiagonalMatrix {
// public:
//     // Member functions

//     // Constructor
//     TridiagonalMatrix(unsigned);

//     double & element(unsigned, unsigned);

//     void fill(double=0.0);

//     void solve(std::vector<double>&);

//     void solve_and_add_result_to_vector(std::vector<double>&, std::vector<double>&);

//     // Member data
//     unsigned dimension;
//     unsigned bandwidth;
//     std::vector<std::vector<double> > data;
// };

// void dxpy(std::vector<double>&, std::vector<double>&);



class DenseMatrix {
public:
    // Member functions

    // Constructor
    DenseMatrix(unsigned);

    double & element(unsigned, unsigned);

    void fill(double=0.0);

    void solve(std::vector<double>&);

    void print();

    double max();

    double min();

    // Member data
    unsigned dimension;
    std::vector<double> data;
};



#endif
