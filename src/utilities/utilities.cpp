#include <vector>
#include <ios>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cmath>

#include "utilities.h"

using std::vector;
using std::streamsize;
using std::setprecision;
using std::cout;
using std::endl;
using std::max_element;
using std::string;
using std::ofstream;

void print_vector(vector<double> vec, unsigned integral_digits,
    unsigned fractional_digits, unsigned linewidth, bool suppress)
{
    typedef vector<double>::size_type vec_sz;

    streamsize prec = cout.precision();

    double max = *std::max_element(vec.begin(), vec.end());
    double min = *std::min_element(vec.begin(), vec.end());

    if (max >= std::pow(10, integral_digits) || min <= -std::pow(10, integral_digits)) {
        cout << std::scientific;
        integral_digits = 5;
    } else {
        cout << std::fixed;
    }

    cout << setprecision(fractional_digits) << '[';

    for (vec_sz i = 0; i != vec.size(); ++i) {

        cout << std::setfill(' ') << std::setw(2 + integral_digits + fractional_digits);

        cout << vec[i];

        if (i != vec.size() - 1) {
           cout << ' ';
        }
    }

    cout << ']' << setprecision(prec) << endl;

    cout << std::defaultfloat;
}

template <class T>
void print_matrix(T mat, unsigned integral_digits, unsigned fractional_digits,
    unsigned linewidth, bool suppress)
{
    streamsize prec = cout.precision();

    auto max = mat.max();
    auto min = mat.min();

    if (max >= std::pow(10.0, integral_digits) || min <= -std::pow(10.0, integral_digits)) {
        cout << std::scientific;

        integral_digits = 5;
    } else {
        cout << std::fixed;
    }

    cout << setprecision(fractional_digits) << '[';

    for (unsigned i = 0; i != mat.dimension; ++i) {

        cout << '[';

        for (unsigned j = 0; j != mat.dimension; ++j) {

            cout << std::setfill(' ') << std::setw(2 + integral_digits + fractional_digits);

            cout << mat.element(i, j);

            if (j != mat.dimension - 1) {
                cout << ' ';
            }
        }

        cout << ']';

        if (i != mat.dimension - 1) {
            cout << endl << ' ';
        }
    }

    cout << ']' << setprecision(prec) << endl;

    cout << std::defaultfloat;
}

void save_vector(std::vector<double> vec, string filename)
{
    ofstream outfile(filename, ofstream::out);

    streamsize prec = outfile.precision();

    outfile << setprecision(10) << std::scientific;

    for (unsigned i = 0; i != vec.size(); i++) {
        outfile << vec[i] << endl;
    }
    outfile << std::fixed << setprecision(prec);

    outfile.close();
}

void save_matrix(std::vector<std::vector<double> > mat, string filename, string delimiter)
{
    ofstream outfile(filename, ofstream::out);

    streamsize prec = outfile.precision();

    outfile << setprecision(10) << std::scientific;

    for (unsigned i = 0; i != mat.size(); i++) {
        for (unsigned j = 0; j != mat[i].size() - 1; j++) {
            outfile << mat[i][j] << delimiter;
        }

        outfile << mat[i].back() << endl;
    }

    outfile << std::fixed << setprecision(prec);

    outfile.close();
}

vector<double> linspace(double start, double stop, unsigned num, bool endpoint)
{
    double h = (stop - start) / static_cast<double>(num - 1);

    vector<double> result(num);

    typename vector<double>::iterator x;

    double val;

    for (x = result.begin(), val = start; x != result.end(); ++x, val += h) {
        *x = val;
    }

    return result;
}
