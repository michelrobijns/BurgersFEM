#ifndef GUARD_tridiagonal_matrix_h
#define GUARD_tridiagonal_matrix_h

#include <vector>
#include <array>

class TridiagonalMatrix {
public:
    // Member functions

    // Constructor
    TridiagonalMatrix(unsigned);

    double & element(unsigned, unsigned);

    void fill(double=0.0);

    void solve(std::vector<double>&);

    //void solve_cyclic(std::vector<double>&, double, double);
    void solve_cyclic(std::vector<double>&);
    
    void print(unsigned=3, unsigned=3, unsigned=100, bool=true);

    // Member data
    unsigned dimension;
    unsigned bandwidth;
    std::vector<std::vector<double> > data;
};

void dxpy(std::vector<double>&, std::vector<double>&);

#endif
