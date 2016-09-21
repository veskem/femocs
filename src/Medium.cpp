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

Medium::Medium() {
    init_statistics();
    reserve(0);
}

// Sort the atoms by their cartesian coordinate or radial coordinate from origin
const void Medium::sort_atoms(const int coord, const string& direction, const Point2 &origin) {
    require(coord >= 0 && coord <= 3, "Invalid coordinate: " + to_string(coord));

    if (coord == 3)
        for(int i = 0; i < get_n_atoms(); ++i)
            atoms[i].point.r = origin.distance2(get_point2(i));

    if (direction == "up" || direction == "asc")
        sort(atoms.begin(), atoms.end(), Atom::sort_up(coord));
    else if (direction == "down" || direction == "desc")
        sort(atoms.begin(), atoms.end(), Atom::sort_down(coord));
}

// Sort the atoms first by first and then by second cartesian coordinate
const void Medium::sort_atoms(const int x1, const int x2, const string& direction, const Point2 &origin) {
    require(x1 >= 0 && x1 <= 3 && x2 >= 0 && x2 <= 3, "Invalid coordinates: " + to_string(x1) + ", " + to_string(x2));

    if (x1 == 3 || x2 == 3)
        for(int i = 0; i < get_n_atoms(); ++i)
            atoms[i].point.r = origin.distance2(get_point2(i));

    if (direction == "up" || direction == "asc")
        sort( atoms.begin(), atoms.end(), Atom::sort_up2(x1, x2) );
    else if (direction == "down" || direction == "desc")
        sort( atoms.begin(), atoms.end(), Atom::sort_down2(x1, x2) );
}

// Reserve memory for data vectors
const void Medium::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid number of atoms: " + n_atoms);
    atoms.reserve(n_atoms);
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
        add_atom(m->get_atom(i));

    this->calc_statistics();
}

// Add atom to atoms vector
const void Medium::add_atom(const Atom& atom) {
    expect(get_n_atoms() <= this->atoms.capacity(), "Allocated vector sizes exceeded!");
    this->atoms.push_back(atom);
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
    return atoms.size();
}

// Get 2-dimensional coordinates of i-th atom
const Point2 Medium::get_point2(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return Point2(atoms[i].point.x, atoms[i].point.y);
}

// Get 3-dimensional coordinates of i-th atom
const Point3 Medium::get_point(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i) + "/" + to_string(get_n_atoms()));
    return atoms[i].point;
}

// Get atom ID
const int Medium::get_id(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return atoms[i].id;
}

// Get atom coordination
const int Medium::get_coordination(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return atoms[i].coord;
}

// Get i-th Atom
const Atom Medium::get_atom(const int i) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return atoms[i];
}

// Set entry to id-s vector
const void Medium::set_id(const int i, const int id) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].id = id;
}

// Set entry to point-s vector
const void Medium::set_point(const int i, const Point3& p) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].point = p;
}

// Set entry to x coordinate vector
const void Medium::set_x(const int i, const double x) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].point.x = x;
}

// Set entry to y coordinate vector
const void Medium::set_y(const int i, const double y) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].point.y = y;
}

// Set entry to z coordinate vector
const void Medium::set_z(const int i, const double z) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].point.z = z;
}

// Set coordination of atom
const void Medium::set_coordination(const int i, const int coord) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    require(coord >= 0, "Invalid coordination: " + to_string(coord));
    atoms[i].coord = coord;
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
    strs << atoms[i];
    return strs.str();
}

} /* namespace femocs */
