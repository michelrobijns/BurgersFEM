#include <iostream>
#include <vector>
#include <algorithm>  // std::fill
#include <Accelerate/Accelerate.h>

#include "linear_algebra.h"

using std::cout;
using std::endl;
using std::vector;
using std::fill;

TridiagonalMatrix::TridiagonalMatrix(unsigned dimension) : dimension(dimension),
    data(vector<vector<double> >(3, vector<double>(dimension, 0.0)))
{
    ;
}

double & TridiagonalMatrix::element(unsigned i, unsigned j)
{
    // Maybe do some bound checking here

    return data[i - j + 1][i];
}

void TridiagonalMatrix::fill(double value)
{
    std::fill(data[0].begin(), data[0].end(), value);
    std::fill(data[1].begin(), data[1].end(), value);
    std::fill(data[2].begin(), data[2].end(), value);
}

void TridiagonalMatrix::solve(vector<double>& d)
{
    vector<double>& a = data[2];  // Subdiagonal
    vector<double>& b = data[1];  // Diagonal
    vector<double>& c = data[0];  // Super diagonal

    unsigned n = dimension - 1;

    c[0] /= b[0];
    d[0] /= b[0];

    for (unsigned i = 1; i != n; ++i) {
        c[i] /= b[i] - a[i] * c[i-1];
        d[i] = (d[i] - a[i] * d[i-1]) / (b[i] - a[i] * c[i-1]);
    }

    d[n] = (d[n] - a[n] * d[n-1]) / (b[n] - a[n] * c[n-1]);

    //for (int i = n; i-- > 0;) {
    for (int i = n - 1; i >= 0; --i) {
        d[i] -= c[i] * d[i+1];
    }

    // int dim = dimension;
    // int nrhs = 1;
    // int info;

    // dgtsv_(&dim, &nrhs, & *data[2].begin()+1, & *data[1].begin(),
    //        & *data[0].begin(), & *d.begin(), &dim, &info);
}

void TridiagonalMatrix::solve_cyclic(vector<double>& d, double alpha, double beta)
{
    vector<double> u(dimension, 0.0);
    vector<double> v(dimension, 0.0);

    double gamma = -data[1][0];

    u.front() = gamma;
    u.back() = alpha;

    v.front() = 1.0;
    v.back() = beta / gamma;

    vector<double> data_copy = data;
    solve(d);
    vector<double> y = d;

    data = data_copy;
    solve(u);
    vector<double> z = u;

    double v_dot_y = v.front() * y.front() + v.back() * y.back()
    double v_dot_z = v.front() * z.front() + v.back() * z.back()

    double factor = -1.0 * v_dot_y / (1.0 + v_dot_z);

    vector<double> x = z;

    cblas_daxpy(x.size(), factor, & *x.begin(), 1.0, & *y.begin(), 1.0);
}

void TridiagonalMatrix::solve_and_add_result_to_vector(vector<double>& d,
    vector<double>& result)
{
    vector<double>& a = data[2];  // Subdiagonal
    vector<double>& b = data[1];  // Diagonal
    vector<double>& c = data[0];  // Super diagonal

    unsigned n = dimension - 1;

    c[0] /= b[0];
    d[0] /= b[0];

    for (unsigned i = 1; i != n; ++i) {
        c[i] /= b[i] - a[i] * c[i-1];
        d[i] = (d[i] - a[i] * d[i-1]) / (b[i] - a[i] * c[i-1]);
    }

    d[n] = (d[n] - a[n] * d[n-1]) / (b[n] - a[n] * c[n-1]);
    result[n] += d[n];

    //for (int i = n; i-- > 0;) {
    for (int i = n - 1; i >= 0; --i) {
        d[i] -= c[i] * d[i+1];
        result[i] += d[i];
    }

    /*int nrhs = 1;
    int info;

    //dgtsv_(&dimension, &nrhs, & *data[2].begin()+1, & *data[1].begin(),
    //       & *data[0].begin(), & *d.begin(), &dimension, &info);

    for (int i = 0; i != d.size(); ++i) {
        result[i] += d[i];
    }*/
}



void dxpy(vector<double>& x, vector<double>& y)
{
    cblas_daxpy(x.size(), 1.0, & *x.begin(), 1.0, & *y.begin(), 1.0);
}



DenseMatrix::DenseMatrix(unsigned dimension) : dimension(dimension),
    data(vector<double>(dimension * dimension))
{
    ;
}

double & DenseMatrix::element(unsigned i, unsigned j)
{
    // Maybe do some bound checking here

    return data[i + j * dimension];  // Column major order
}

void DenseMatrix::fill(double value)
{
    std::fill(data.begin(), data.end(), value);
}

void DenseMatrix::solve(vector<double>& d)
{
    // int dim = dimension;
    // int nrhs = 1;
    // int info;

    // dgtsv_(&dim, &nrhs, & *data[2].begin()+1, & *data[1].begin(),
    //        & *data[0].begin(), & *d.begin(), &dim, &info);

    int dim = dimension;
    int nrhs = 1;
    int info;
    vector<int> ipiv(dimension);

    dgesv_(&dim,
           &nrhs,
           & *data.begin(),
           &dim,
           & *ipiv.begin(),
           & *d.begin(),
           &dim,
           &info);

    // int dgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_doublereal *__a,
    //         __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_doublereal *__b,
    //         __CLPK_integer *__ldb,
    //         __CLPK_integer *__info)
}

void DenseMatrix::print()
{
    cout << "Hello, World!" << endl;
}

double DenseMatrix::max()
{
    return *std::max_element(data.begin(), data.end());
}

double DenseMatrix::min()
{
    return *std::min_element(data.begin(), data.end());
}
