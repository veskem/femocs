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
    require(n_atoms >= 0, "Invalid number of atoms: " + n_atoms);
    this->id.reserve(n_atoms);
    this->point.reserve(n_atoms);
    this->coordination.reserve(n_atoms);
}

// Define the addition of two Mediums
Medium& Medium::operator +=(Medium &m) {
    add(&m);
    return *this;
}

// Add data from other Medium to current one
const void Medium::add(Medium *m) {
    const int n_atoms1 = get_n_atoms();
    const int n_atoms2 = m->get_n_atoms();

    this->reserve(n_atoms1 + n_atoms2);

    for(int i = 0; i < n_atoms2; ++i)
        add_atom(m->get_id(i), m->get_point(i), m->get_coordination(i));

    this->calc_statistics();
}

const void Medium::add_atom(const int id, const Point3 &point, const int coord) {
    expect(get_n_atoms() <= this->point.capacity(), "Allocated vector sizes exceeded!");
    this->id.push_back(id);
    this->point.push_back(point);
    this->coordination.push_back(coord);
}

// Initialize statistics about Medium
const void Medium::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
    sizes.xmean = sizes.ymean = sizes.zmean = 0.0;

    sizes.xbox = sizes.ybox = sizes.zbox = 0;
    sizes.zminbox = DBL_MAX;
    sizes.zmaxbox = DBL_MIN;
}

// Calculate the statistics about Medium
const void Medium::calc_statistics() {
    double xx, yy, zz;
    int n_atoms = get_n_atoms();
    init_statistics();

    Point3 average(0,0,0);

    // Find min and max coordinates
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = get_point(i);
        average += point;
        sizes.xmax = max(sizes.xmax, point.x);
        sizes.xmin = min(sizes.xmin, point.x);
        sizes.ymax = max(sizes.ymax, point.y);
        sizes.ymin = min(sizes.ymin, point.y);
        sizes.zmax = max(sizes.zmax, point.z);
        sizes.zmin = min(sizes.zmin, point.z);
    }

    sizes.xmean = average.x / n_atoms;
    sizes.ymean = average.y / n_atoms;
    sizes.zmean = average.z / n_atoms;

    // Define size of simubox
    sizes.xbox = sizes.xmax - sizes.xmin;
    sizes.ybox = sizes.ymax - sizes.ymin;
    sizes.zbox = sizes.zmax - sizes.zmin;
    sizes.zminbox = sizes.zmin;
    sizes.zmaxbox = sizes.zmax;
}

// Get number of atoms in Medium
const int Medium::get_n_atoms() {
    return point.size();
}

// Get 2-dimensional coordinates of i-th atom
const Point2 Medium::get_point2(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    return Point2(point[i].x, point[i].y);
}

// Get 3-dimensional coordinates of i-th atom
const Point3 Medium::get_point(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    return point[i];
}

// Get atom ID
const int Medium::get_id(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    return id[i];
}

// Get atom coordination
const int Medium::get_coordination(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    return coordination[i];
}

// Set entry to x coordinate vector
const void Medium::set_x(const int i, const double x) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    point[i].x = x;
}

// Set entry to y coordinate vector
const void Medium::set_y(const int i, const double y) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    point[i].y = y;
}

// Set entry to z coordinate vector
const void Medium::set_z(const int i, const double z) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    point[i].z = z;
}

// Set coordination of atom
const void Medium::set_coordination(const int i, const int coord) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + i);
    require(coord >= 0, "Invalid coordination: " + coord);
    this->coordination[i] = coord;
}

// Pick the suitable write function based on the file type
// Function works only in debug mode
void Medium::output(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string ftype = get_file_type(file_name);
    expect(ftype == "xyz", "Unsupported file type: " + ftype);

    const int n_atoms = get_n_atoms();

    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << n_atoms << "\n";
    out_file << get_data_string(-1) << endl;

    for (int i = 0; i < n_atoms; ++i)
        out_file << get_data_string(i) << endl;

    out_file.close();
}

// Compile data string from the data vectors
const string Medium::get_data_string(const int i) {
    if(i < 0) return "Medium data: id x y z coordination";

    ostringstream strs;
    strs << get_id(i) << " " << get_point(i) << " " << get_coordination(i);
    return strs.str();
}

} /* namespace femocs */
