/*
 * Medium.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Surface.h"
#include <float.h>

#include <fstream>

using namespace std;

namespace femocs {

// Constructor initializes arrays
Surface::Surface(const int natoms, double latconst) {
    reserve(natoms);
    init_statistics();
    sizes.latconst = latconst;
}
;

// Add atom coordinates
const void Surface::add_atom(const double x, const double y, const double z, const int coord,
        const int type) {
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->coordination.push_back(coord);
    this->type.push_back(type);
}

const void Surface::reserve(const int N) {
    this->x.reserve(N);                // Default values are 0
    this->y.reserve(N);
    this->z.reserve(N);
    this->Sx.reserve(N);
    this->Sy.reserve(N);
    this->Sz.reserve(N);
    this->coordination.reserve(N);
    this->type.reserve(N);
    this->isEvaporated.reserve(N);    // Default values are False
}

vector<double>* Surface::get_xs() {
    return &x;
}
vector<double>* Surface::get_ys() {
    return &y;
}
vector<double>* Surface::get_zs() {
    return &z;
}

void Surface::set_x(const int i, const double x) {
    this->x[i] = x;
}

void Surface::set_y(const int i, const double y) {
    this->y[i] = y;
}

void Surface::set_z(const int i, const double z) {
    this->z[i] = z;
}

// Get element from x, y or z vector
const double Surface::get_x(const int i) {
    return x[i];
}
const double Surface::get_y(const int i) {
    return y[i];
}
const double Surface::get_z(const int i) {
    return z[i];
}

// Get atom type
const int Surface::get_type(const int i) {
    return type[i];
}
// Get atom coordination
const int Surface::get_coordination(const int i) {
    return coordination[i];
}

// Get number of atoms on surface
const int Surface::get_N() {
    return x.size();
}

const void Surface::calc_statistics() {
    init_statistics();
    for (int i = 0; i < get_N(); ++i) {
        if(sizes.xmax < get_x(i)) sizes.xmax = get_x(i);
        if(sizes.xmin > get_x(i)) sizes.xmin = get_x(i);
        if(sizes.ymax < get_y(i)) sizes.ymax = get_y(i);
        if(sizes.ymin > get_y(i)) sizes.ymin = get_y(i);
        if(sizes.zmax < get_z(i)) sizes.zmax = get_z(i);
        if(sizes.zmin > get_z(i)) sizes.zmin = get_z(i);
    }
}

const void Surface::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
}

// Output surface data to file
const void Surface::output(const string& file_name) {
    int nrOfAtoms = this->x.size();
    ofstream myfile;
    myfile.open(file_name);
    myfile << nrOfAtoms << "\n";
    myfile << "Surface of the nanotip\n";

    for (int i = 0; i < nrOfAtoms; ++i) {
        myfile << i << " ";
        myfile << this->x[i] << " ";
        myfile << this->y[i] << " ";
        myfile << this->z[i] << " ";
        myfile << "1" << " ";
        myfile << this->coordination[i] << endl;
    }
    myfile.close();
}

} /* namespace femocs */
