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
    this->id.reserve(n_atoms);
    this->point.reserve(n_atoms);
    this->coordination.reserve(n_atoms);
}

const void Medium::add_atom(const int id, const Point3d &point, const int coord) {
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
}

// Calculate the statistics about Medium
const void Medium::calc_statistics() {
    double xx, yy, zz;
    int n_atoms = get_n_atoms();
    init_statistics();

    Point3d average(0,0,0);

    // Find min and max coordinates
    for (int i = 0; i < n_atoms; ++i) {
        Point3d point = get_point(i);
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
}

// Get number of atoms in Medium
int Medium::get_n_atoms() {
    return point.size();
}

// Get i-th 2-dimensional point
const Point2d Medium::get_point2d(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return Point2d(point[i].x, point[i].y);
}

// Get i-th 3-dimensional point
Point3d Medium::get_point(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return point[i];
}

// Get entry from ID-s vector
const int Medium::get_id(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return id[i];
}

// Get atom coordination
const int Medium::get_coordination(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return coordination[i];
}

// Set entry to x coordinate vector
void Medium::set_x(const int i, const double x) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    point[i].x = x;
}

// Set entry to y coordinate vector
void Medium::set_y(const int i, const double y) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    point[i].y = y;
}

// Set entry to z coordinate vector
void Medium::set_z(const int i, const double z) {
    require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    point[i].z = z;
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
    strs << get_id(i) << " " << get_point(i) << " " << get_coordination(i);
    return strs.str();
}

} /* namespace femocs */
