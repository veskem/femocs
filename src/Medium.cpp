/*
 * Medium.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Medium.h"

#include <fstream>
#include <sstream>
#include <float.h>

using namespace std;
namespace femocs {

Medium::Medium(){
    init_statistics();
    reserve(0);
}

// Reserve memory for data vectors
const void Medium::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid number of atoms!");
    this->x.reserve(n_atoms);
    this->y.reserve(n_atoms);
    this->z.reserve(n_atoms);
    this->coordination.reserve(n_atoms);
}

// Add atom's parameters to data vectors
const void Medium::add_atom(const double x, const double y, const double z, const int coord) {
    expect(get_n_atoms() <= this->x.capacity(), "Allocated vector sizes exceeded!");
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->coordination.push_back(coord);
}

// Initialize statistics about Medium
const void Medium::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
}

// Calculate the statistics about Medium
const void Medium::calc_statistics() {
    init_statistics();

    // Find min and max coordinates
    for (int i = 0; i < get_n_atoms(); ++i) {
        if (sizes.xmax < get_x(i)) sizes.xmax = get_x(i);
        if (sizes.xmin > get_x(i)) sizes.xmin = get_x(i);
        if (sizes.ymax < get_y(i)) sizes.ymax = get_y(i);
        if (sizes.ymin > get_y(i)) sizes.ymin = get_y(i);
        if (sizes.zmax < get_z(i)) sizes.zmax = get_z(i);
        if (sizes.zmin > get_z(i)) sizes.zmin = get_z(i);
    }
}

// Get number of atoms in Medium
int Medium::get_n_atoms() {
    return x.size();
}

// Get i-th 2-dimensional point
const Point2d Medium::get_point2d(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return Point2d(x[i], y[i]);
}

// Get i-th 3-dimensional point
Point3d Medium::get_point(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return Point3d(x[i], y[i], z[i]);
}

// Get entry from x coordinate vector
const double Medium::get_x(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return x[i];
}

// Get entry from y coordinate vector
const double Medium::get_y(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return y[i];
}

// Get entry from z coordinate vector
const double Medium::get_z(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return z[i];
}

// Get atom coordination
const int Medium::get_coordination(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return coordination[i];
}

// Set entry to x coordinate vector
void Medium::set_x(const int i, const double x) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    this->x[i] = x;
}

// Set entry to y coordinate vector
void Medium::set_y(const int i, const double y) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    this->y[i] = y;
}

// Set entry to z coordinate vector
void Medium::set_z(const int i, const double z) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    this->z[i] = z;
}

// Set coordination of atom
void Medium::set_coordination(const int i, const int coord) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    require(coord >= 0, "Invalid coordination!");
    this->coordination[i] = coord;
}

// Pick the suitable write function based on the file type
// Function works only in debug mode
const void Medium::output(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string ftype = get_file_type(file_name);
    expect(ftype == "xyz", "Unsupported file type!");

    int n_atoms = get_n_atoms();

    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << n_atoms << "\n";
    out_file << get_data_string(-1) << endl;

    for (int i = 0; i < n_atoms; ++i)
        out_file << get_data_string(i) << endl;

    out_file.close();
}

// Extract file extension from file name
const string Medium::get_file_type(const string file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    return file_name.substr(start, end);
}

// Compile data string from the data vectors
const string Medium::get_data_string(const int i) {
    if(i < 0)
        return "Types of data: id x y z coordination";

    ostringstream strs;
    strs << i << " " << get_point(i) << " " << get_coordination(i);
    return strs.str();
}

} /* namespace femocs */
