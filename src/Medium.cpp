/*
 * Medium.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Medium.h"
#include <CGAL/hilbert_sort.h>
#include <float.h>
#include <fstream>
#include <numeric>

using namespace std;
namespace femocs {

Medium::Medium() {
    init_statistics();
    reserve(0);
}

Medium::Medium(const int n_atoms) {
    init_statistics();
    reserve(n_atoms);
}

// Sort the atoms by their cartesian or radial coordinate
const void Medium::sort_atoms(const int coord, const string& direction) {
    require(coord >= 0 && coord <= 3, "Invalid coordinate: " + to_string(coord));

    if (coord == 3) {
        calc_statistics();
        Point2 origin(sizes.xmid, sizes.ymid);
        for (int i = 0; i < get_n_atoms(); ++i)
            set_marker(i, 10000 * origin.distance2(get_point2(i)));
        if (direction == "up" || direction == "asc")
            sort(atoms.begin(), atoms.end(), Atom::sort_marker_up());
        else
            sort(atoms.begin(), atoms.end(), Atom::sort_marker_down());
    }

    else if (direction == "up" || direction == "asc")
        sort(atoms.begin(), atoms.end(), Atom::sort_up(coord));
    else if (direction == "down" || direction == "desc")
        sort(atoms.begin(), atoms.end(), Atom::sort_down(coord));
}

// Sort the atoms first by 1st and then by 2nd cartesian coordinate
const void Medium::sort_atoms(const int x1, const int x2, const string& direction) {
    require(x1 >= 0 && x1 <= 2 && x2 >= 0 && x2 <= 2, "Invalid coordinates: " + to_string(x1) + ", " + to_string(x2));

    if (direction == "up" || direction == "asc")
        sort( atoms.begin(), atoms.end(), Atom::sort_up2(x1, x2) );
    else if (direction == "down" || direction == "desc")
        sort( atoms.begin(), atoms.end(), Atom::sort_down2(x1, x2) );
}

// Perform spatial sorting by ordering atoms along Hilbert curve
const void Medium::sort_spatial() {
    CGAL::hilbert_sort( atoms.begin(), atoms.end(), Atom::sort_spatial(), CGAL::Hilbert_sort_middle_policy() );
}

// Reserve memory for data vectors
const void Medium::reserve(const int n_atoms) {
    require(n_atoms >= 0, "Invalid number of atoms: " + to_string(n_atoms));
    atoms.clear();
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

    atoms.reserve(n_atoms1 + n_atoms2);
    for(int i = 0; i < n_atoms2; ++i)
        add_atom(m->get_atom(i));

    this->calc_statistics();
}

// Add atom to atoms vector
const void Medium::add_atom(const Atom& atom) {
    expect(get_n_atoms() < atoms.capacity(), "Allocated vector sizes exceeded!");
    atoms.push_back(atom);
}

// Add atom with defalt id and marker to atoms vector
const void Medium::add_atom(const Point3& point) {
    expect(get_n_atoms() < atoms.capacity(), "Allocated vector sizes exceeded!");
    atoms.push_back(Atom(-1, point, 0));
}

// Initialize statistics about Medium
const void Medium::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax =-DBL_MAX;
    sizes.xmean = sizes.ymean = sizes.zmean = 0.0;
    sizes.xmid = sizes.ymid = sizes.zmid = 0.0;

    sizes.xbox = sizes.ybox = sizes.zbox = 0;
    sizes.zminbox = DBL_MAX;
    sizes.zmaxbox =-DBL_MAX;
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
        sizes.xmin = min(sizes.xmin, point.x);
        sizes.xmax = max(sizes.xmax, point.x);
        sizes.ymin = min(sizes.ymin, point.y);
        sizes.ymax = max(sizes.ymax, point.y);
        sizes.zmin = min(sizes.zmin, point.z);
        sizes.zmax = max(sizes.zmax, point.z);
    }

    // Define average coordinates
    average *= (1.0 / n_atoms);
    sizes.xmean = average.x;
    sizes.ymean = average.y;
    sizes.zmean = average.z;

    // Define size of simubox
    sizes.xbox = sizes.xmax - sizes.xmin;
    sizes.ybox = sizes.ymax - sizes.ymin;
    sizes.zbox = sizes.zmax - sizes.zmin;
    sizes.zminbox = sizes.zmin;
    sizes.zmaxbox = sizes.zmax;

    // Define the centre of simubox
    sizes.xmid = (sizes.xmax + sizes.xmin) / 2;
    sizes.ymid = (sizes.ymax + sizes.ymin) / 2;
    sizes.zmid = (sizes.zmax + sizes.zmin) / 2;
}

// Get number of atoms in Medium
const int Medium::get_n_atoms() const {
    return atoms.size();
}

// Get 2-dimensional coordinates of i-th atom
const Point2 Medium::get_point2(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return Point2(atoms[i].point.x, atoms[i].point.y);
}

// Get 3-dimensional coordinates of i-th atom
const Point3 Medium::get_point(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i) + "/" + to_string(get_n_atoms()));
    return atoms[i].point;
}

// Get atom ID
const int Medium::get_id(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return atoms[i].id;
}

// Get atom marker
const int Medium::get_marker(const int i) const {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return atoms[i].marker;
}

// Get i-th Atom
const Atom Medium::get_atom(const int i) const {
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

// Set atom marker
const void Medium::set_marker(const int i, const int m) {
    require(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    atoms[i].marker = m;
}

// Pick the suitable write function based on the file type
// Function works only in debug mode
void Medium::write(const string &file_name, const int n_max) const {
    if (!MODES.WRITEFILE) return;
    
    int n_atoms = get_n_atoms();
    if (n_max > 0 && n_max < n_atoms) n_atoms = n_max;
    
    expect(n_atoms > 0, "Zero atoms detected!");
    string ftype = get_file_type(file_name);
    
    ofstream outfile;
    if (ftype == "movie") outfile.open(file_name, ios_base::app);
    else outfile.open(file_name);
    require(outfile.is_open(), "Can't open a file " + file_name);
    
    if (ftype == "xyz" || ftype == "movie")
        write_xyz(outfile, n_atoms);
    else if (ftype == "vtk")
        write_vtk(outfile, n_atoms);
    else    
        require(false, "Unsupported file type: " + ftype);

    outfile.close();
}

// Compile data string from the data vectors
const string Medium::get_data_string(const int i) const {
    if(i < 0) return "Medium properties=id:R:1:pos:R:3:marker:R:1";

    ostringstream strs;
    strs << atoms[i];
    return strs.str();
}

// Output atom data in .xyz format
const void Medium::write_xyz(ofstream& out, const int n_atoms) const {
    out << n_atoms << "\n";
    out << get_data_string(-1) << endl;

    for (int i = 0; i < n_atoms; ++i)
        out << get_data_string(i) << endl;
}

// Output atom data in .vtk format
const void Medium::write_vtk(ofstream& out, const int n_atoms) const {
    out << "# vtk DataFile Version 3.0\n";
    out << "# Medium data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_atoms << " double\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << get_point(i) << "\n";

    get_cell_types(out, n_atoms);

    out << "\nPOINT_DATA " << n_atoms << "\n";
    get_point_data(out, n_atoms);

    out << "\nCELL_DATA " << n_atoms << "\n";
    get_cell_data(out, n_atoms);
}

// Get point representation in vtk format
const void Medium::get_cell_types(ofstream& out, const int n_cells) const {
    const int dim = 1;
    const int celltype = 1; // cell == vertex

    // Output the vertices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (size_t i = 0; i < n_cells; ++i)
        out << "1 " << i << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (size_t i = 0; i < n_cells; ++i)
        out << celltype << "\n";
}

// Get data scalar and vector data associated with vtk cells
const void Medium::get_cell_data(ofstream& out, const int n_cells) const {}

// Get data scalar and vector data associated with vtk nodes
const void Medium::get_point_data(ofstream& out, const int n_atoms) const {
    // write IDs of atoms
    out << "SCALARS id int\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << atoms[i].id << "\n";

    // write atom markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << "\n";
}

} /* namespace femocs */
