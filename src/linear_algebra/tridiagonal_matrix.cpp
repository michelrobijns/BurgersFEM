#include <iostream>  // std::cout, std::endl;
#include <vector>  // std::vector
#include <algorithm>  // std::fill
#include <Accelerate/Accelerate.h>  // BLAS and LAPACK routines

#include <ios>
#include <iomanip>
#include <iostream>
#include <string>

#include "tridiagonal_matrix.h"

using std::cout;
using std::endl;
using std::vector;
using std::fill;

using std::streamsize;
using std::setprecision;
using std::max_element;
using std::string;


TridiagonalMatrix::TridiagonalMatrix(unsigned dimension) : dimension(dimension),
    data(vector<vector<double> >(3, vector<double>(dimension, 0.0))) {}

double & TridiagonalMatrix::element(unsigned i, unsigned j)
{
    // Maybe do some bound checking here

    return data[i - j + 1][i];
}

void TridiagonalMatrix::fill(double value)
{
    std::fill(data[0].begin(), data[0].end()-1, value);
    std::fill(data[1].begin(), data[1].end(), value);
    std::fill(data[2].begin()+1, data[2].end(), value);
}

void TridiagonalMatrix::solve(vector<double>& b)
{
    // Thomas algorithm

    vector<double>& sub = data[2];    // Subdiagonal
    vector<double>& diag = data[1];   // Diagonal
    vector<double>& super = data[0];  // Superdiagonal

    unsigned n = dimension - 1;

    super[0] /= diag[0];
    b[0] /= diag[0];

    for (unsigned i = 1; i != n; ++i) {
        super[i] /= diag[i] - sub[i] * super[i-1];
        b[i] = (b[i] - sub[i] * b[i-1]) / (diag[i] - sub[i] * super[i-1]);
    }

    b[n] = (b[n] - sub[n] * b[n-1]) / (diag[n] - sub[n] * super[n-1]);

    for (int i = n - 1; i >= 0; --i) {
        b[i] -= super[i] * b[i+1];
    }

    // int dim = dimension;
    // int nrhs = 1;
    // int info;

    // dgtsv_(&dim, &nrhs, & *data[2].begin()+1, & *data[1].begin(),
    //        & *data[0].begin(), & *d.begin(), &dim, &info);
}

void TridiagonalMatrix::solve_cyclic(vector<double>& b)
{
    element(0, 0) += element(dimension-1, dimension-1);
    b[0] += b[dimension-1];

    // Set `alpha', `beta', and `gamma'
    double alpha = element(dimension-2, dimension-1);
    double beta = element(dimension-1, dimension-2);
    double gamma = -1.0 * element(0, 0);

    // "remove" the last row and the last column of the matrix
    element(dimension-1, dimension-1) = 1.0;
    element(dimension-2, dimension-1) = 0.0;
    element(dimension-1, dimension-2) = 0.0;

    // "remove" the right-hand vector
    b.back() = 0.0;

    // This is where the Sherman-Morrison formula starts
    element(0, 0) -= gamma;
    element(dimension-2, dimension-2) -= alpha * beta / gamma;

    // Allocate additional storage
    vector<double> u(dimension, 0.0);
    vector<double> v(dimension, 0.0);
    vector<vector<double> > data_copy = data;

    // Set `u'
    u[0] = gamma;
    u[dimension-2] = alpha;

    // Set 'v'
    v[0] = 1.0;
    v[dimension-2] = beta / gamma;

    // Solve A * y = b and store the result in `b'
    solve(b);

    // Solve A * z = u and store the result in `u'
    data = data_copy;
    solve(u);

    double v_dot_y = cblas_ddot(dimension, & *v.begin(), 1.0, & *b.begin(), 1.0);
    double v_dot_z = cblas_ddot(dimension, & *v.begin(), 1.0, & *u.begin(), 1.0);

    double factor = -1.0 * v_dot_y / (1.0 + v_dot_z);

    // Compute x = y + factor * b and store the result in `b'
    cblas_daxpy(dimension, factor, & *u.begin(), 1.0, & *b.begin(), 1.0);

    b[dimension-1] = b[0];
}

// void TridiagonalMatrix::solve_cyclic(vector<double>& b, double alpha, double beta)
// {
//     double gamma = -1.0 * element(0, 0);

//     // Modify matrix
//     element(0, 0) -= gamma;
//     element(dimension-1, dimension-1) -= alpha * beta / gamma;

//     // Allocate additional storage
//     vector<double> u(dimension, 0.0);
//     vector<double> v(dimension, 0.0);
//     vector<vector<double> > data_copy = data;

//     // Set `u'
//     u.front() = gamma;
//     u.back() = alpha;

//     // Set 'v'
//     v.front() = 1.0;
//     v.back() = beta / gamma;

//     // Solve A * y = b and store the result in `b'
//     solve(b);

//     // Solve A * z = u and store the result in `u'
//     data = data_copy;
//     solve(u);

//     double v_dot_y = cblas_ddot(dimension, & *v.begin(), 1.0, & *b.begin(), 1.0);
//     double v_dot_z = cblas_ddot(dimension, & *v.begin(), 1.0, & *u.begin(), 1.0);

//     double factor = -1.0 * v_dot_y / (1.0 + v_dot_z);

//     // Compute x = y + factor * b and store the result in `b'
//     cblas_daxpy(dimension, factor, & *u.begin(), 1.0, & *b.begin(), 1.0);
// }

void TridiagonalMatrix::print(unsigned integral_digits, unsigned fractional_digits, unsigned linewidth, bool suppress)
{
    streamsize prec = cout.precision();

    auto max = 3.0;   // Change
    auto min = -1.0;  // change

    if (max >= std::pow(10.0, integral_digits) || min <= -std::pow(10.0, integral_digits)) {
        cout << std::scientific;

        integral_digits = 5;
    } else {
        cout << std::fixed;
    }

    cout << setprecision(fractional_digits) << '[';

    for (unsigned i = 0; i != dimension; ++i) {

        cout << '[';

        for (unsigned j = 0; j != dimension; ++j) {

            cout << std::setfill(' ') << std::setw(2 + integral_digits + fractional_digits);

            if (i == j || j == (i + 1) || j == (i - 1)) {
                cout << element(i, j);
            } else {
                cout << 0.0;
            }

            if (j != dimension - 1) {
                cout << ' ';
            }
        }

        cout << ']';

        if (i != dimension - 1) {
            cout << endl << ' ';
        }
    }

    cout << ']' << setprecision(prec) << endl;

    cout << std::defaultfloat;
}

void dxpy(vector<double>& x, vector<double>& y)
{
    cblas_daxpy(x.size(), 1.0, & *x.begin(), 1.0, & *y.begin(), 1.0);
}
