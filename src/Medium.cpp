/*
 * Medium.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Medium.h"

#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;
namespace femocs {

// Function to calculate the indices of sorted array
template<class Vals>
void get_sort_permutation(const Vals& values, vector<int>* v) {
    int size = values.size();
    v->clear();
    v->reserve(size);
    for (int i = 0; i < size; ++i)
        v->push_back(i);

    sort(v->begin(), v->end(), [&values](int a, int b) -> bool {
        return values[a] < values[b];
    });
}

// Get element from x, y or z vector
const double Medium::get_x(const int i) {
    return x[i];
}
const double Medium::get_y(const int i) {
    return y[i];
}
const double Medium::get_z(const int i) {
    return z[i];
}
// Get atom type
const int Medium::get_type(const int i) {
    return type[i];
}
// Get number of atoms on surface
const int Medium::get_n_atoms() {
    return N;
}
// Output surface data to file
const void Medium::output(const string& file_name) {
    int nrOfAtoms = this->x.size();
    ofstream myfile;
    myfile.open(file_name);
    myfile << nrOfAtoms << "\n";
    myfile << "Data of the nanotip: id x y z type\n";

    for (int i = 0; i < nrOfAtoms; ++i) {
        myfile << i << " ";
        myfile << this->x[i] << " ";
        myfile << this->y[i] << " ";
        myfile << this->z[i] << " ";
        myfile << this->type[i] << endl;
    }
    myfile.close();
}
// Function to sort the surface atoms according to their radial distance on xy-plane
void Medium::sort_atoms(vector<int>* permutation_indxs) {
    int i;
    double xloc, yloc;
    vector<double> r2;
    r2.reserve(N);

    for (i = 0; i < N; ++i) {
        xloc = this->x[i];
        yloc = this->y[i];
        r2.push_back(xloc * xloc + yloc * yloc);
    }
    get_sort_permutation(r2, permutation_indxs);

    for (i = 0; i < N; ++i)
        cout << r2[i] << endl;
}

} /* namespace femocs */

