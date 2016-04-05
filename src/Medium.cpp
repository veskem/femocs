/*
 * Medium.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Medium.h"

using namespace std;
namespace femocs {

// Initialize statistics about Medium
const void Medium::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
}

// Calculate the statistics about Medium
const void Medium::calc_statistics() {
    init_statistics();
    for (int i = 0; i < get_n_atoms(); ++i) {
        if(sizes.xmax < get_x(i)) sizes.xmax = get_x(i);
        if(sizes.xmin > get_x(i)) sizes.xmin = get_x(i);
        if(sizes.ymax < get_y(i)) sizes.ymax = get_y(i);
        if(sizes.ymin > get_y(i)) sizes.ymin = get_y(i);
        if(sizes.zmax < get_z(i)) sizes.zmax = get_z(i);
        if(sizes.zmin > get_z(i)) sizes.zmin = get_z(i);
    }
}

// Get entry from x coordinate vector
const double Medium::get_x(const int i) {
    return x[i];
}

// Get entry from y coordinate vector
const double Medium::get_y(const int i) {
    return y[i];
}

// Get entry from z coordinate vector
const double Medium::get_z(const int i) {
    return z[i];
}

// Get atom coordination
const int Medium::get_coordination(const int i) {
    return coordination[i];
}

// Get number of atoms in Medium
const int Medium::get_n_atoms() {
    return x.size();
}

// Set entry to x coordinate vector
void Medium::set_x(const int i, const double x) {
    this->x[i] = x;
}

// Set entry to y coordinate vector
void Medium::set_y(const int i, const double y) {
    this->y[i] = y;
}

// Set entry to z coordinate vector
void Medium::set_z(const int i, const double z) {
    this->z[i] = z;
}

// Set coordination of atom
void Medium::set_coordination(const int i, const int coord) {
    this->coordination[i] = coord;
}

// Output surface data to file in xyz format
const void Medium::output(const string& file_name) {
    int n_atoms = this->x.size();
    ofstream myfile;
    myfile.open(file_name);
    myfile << n_atoms << "\n";
    myfile << "Data of the nanotip: id x y z type\n";

    for (int i = 0; i < n_atoms; ++i)
        myfile << get_data_string(i) << endl;
    
    myfile.close();
}

// Compile data string from the data vectors
const string Medium::get_data_string(const int i) {
    ostringstream strs;
    strs << i << " " << get_x(i) << " " << get_y(i) << " " << get_z(i) <<" "<< get_coordination(i);
    return strs.str();
}

const void Medium::reserve(const int n_atoms) {
    this->x.reserve(n_atoms);
    this->y.reserve(n_atoms);
    this->z.reserve(n_atoms);
    this->coordination.reserve(n_atoms);
}

const void Medium::add_atom(const double x, const double y, const double z, const int coord) {
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->coordination.push_back(coord);
}

} /* namespace femocs */

