/*
 * AtomReader.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "AtomReader.h"

#include <stddef.h>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
namespace femocs {

// AtomReader constructor
AtomReader::AtomReader() {
    init_statistics();
}

// Initialise statistics about coordinates
const void AtomReader::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
}

// Calculate the statistics about coordinates
const void AtomReader::calc_statistics() {
    init_statistics();

    // Find min and max of coordinates
    for (int i = 0; i < get_n_atoms(); ++i) {
        if (sizes.xmax < get_x(i)) sizes.xmax = get_x(i);
        if (sizes.xmin > get_x(i)) sizes.xmin = get_x(i);
        if (sizes.ymax < get_y(i)) sizes.ymax = get_y(i);
        if (sizes.ymin > get_y(i)) sizes.ymin = get_y(i);
        if (sizes.zmax < get_z(i)) sizes.zmax = get_z(i);
        if (sizes.zmin > get_z(i)) sizes.zmin = get_z(i);
    }

    // Define size of simubox
    sizes.xbox = sizes.xmax - sizes.xmin;
    sizes.ybox = sizes.ymax - sizes.ymin;
    sizes.zbox = sizes.zmax - sizes.zmin;
    sizes.zminbox = sizes.zmin;
    sizes.zmaxbox = sizes.zmax;
}

// Reserve memory for data vectors
const void AtomReader::reserve(const int n_atoms) {
	require(n_atoms > 0, "Invalid number of atoms!");
    x.reserve(n_atoms);
    y.reserve(n_atoms);
    z.reserve(n_atoms);
    type.reserve(n_atoms);
}

// Add atom with its coordinates and type
const void AtomReader::add_atom(const double x, const double y, const double z, const int type) {
	expect(get_n_atoms() < this->x.capacity(), "Allocated vector sizes exceeded!");
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->type.push_back(type);
}

// Calculate coordination (nnn within cutoff radius) of all the atoms
// Pick suitable algorithm depending simulation type (MD or KMC)
const void AtomReader::calc_coordination(const double cutoff, const int nnn) {
	require(cutoff > 0 && nnn >= 0, "Invalid cutoff or nnn!");
    if (types.simu_type == "md")
        calc_md_coordination(cutoff, nnn);
    else if (types.simu_type == "kmc")
    	calc_kmc_coordination(nnn);
}

// Calculate coordination for atoms coming from MD simulations
const void AtomReader::calc_md_coordination(const double cutoff, const int nnn) {
	require(cutoff > 0 && nnn >= 0, "Invalid cutoff or nnn!");
    int N = get_n_atoms();
    int i, j, coord;
    double r2, dx, dy, dz;
    double cutoff2 = cutoff * cutoff;

    coordination.reserve(N);

//    omp_set_num_threads(this->nrOfOmpThreads);
//#pragma omp parallel for shared(coords,is_surface) private(i,j,dx,dy,dz,r2) reduction(+:coord,N)
    for (i = 0; i < N; ++i) {
        coord = 0;
        //!< TODO Implement neighborlist here!
        for (j = 0; j < N; ++j) {
            if (i == j) continue;
            dx = fabs(x[i] - x[j]);
            dy = fabs(y[i] - y[j]);
            dz = fabs(z[i] - z[j]);

            dx = min(dx, fabs(dx - sizes.xbox)); // apply periodic boundary condition in x-direction
            dy = min(dy, fabs(dy - sizes.ybox)); // apply periodic boundary condition in y-direction
            //dz = min(dz, fabs(dz - sizes.zbox)); // apply periodic boundary condition in z-direction

            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 <= cutoff2) coord++;
            if (coord >= nnn) break; // Coordination can't be bigger than the biggest expected
        }

        coordination[i] = coord;
    }
}

// Calculate coordination for atoms coming from KMC simulations
const void AtomReader::calc_kmc_coordination(const int nnn) {
	require(nnn >= 0, "Invalid nnn!");
    int N = get_n_atoms();
    int i, j, coord;

    coordination.reserve(N);

    for (i = 0; i < N; ++i) {
        if (type[i] == types.type_bulk)
            coordination[i] = nnn;
        else if (type[i] == types.type_surf)
            coordination[i] = (int) nnn / 2;
        else
            coordination[i] = 0;
    }
}

// Redefine simubox size for example to insert electric field height
const void AtomReader::resize_box(const double zmin, const double zmax) {
	require(zmin <= zmax, "Invalid size for simulation box!");
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;
    sizes.zbox = zmax - zmin;
}

// =================================
// *** GETTERS: ***************
const int AtomReader::get_n_atoms() {
    return x.size();
}

const double AtomReader::get_x(const int i) {
	require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return x[i];
}

const double AtomReader::get_y(const int i) {
	require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return y[i];
}

const double AtomReader::get_z(const int i) {
	require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return z[i];
}

const int AtomReader::get_type(const int i) {
	require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return type[i];
}

const int AtomReader::get_coord(const int i) {
	require(i >= 0 && i < get_n_atoms(), "Invalid index!");
    return coordination[i];
}

const string AtomReader::get_file_type(const string file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    string file_type = file_name.substr(start, end);

    if (file_type == "xyz") this->types.simu_type = "md";
    if (file_type == "dump") this->types.simu_type = "md";
    if (file_type == "ckx") this->types.simu_type = "kmc";

    return file_type;
}

// =================================
// *** IMPORTERS: ***************

const void AtomReader::import_helmod() {
	require(false, "AtomReader::import_helmod() not implemented!");
}

const void AtomReader::import_kimocs() {
	require(false, "AtomReader::import_kimocs() not implemented!");
}

const void AtomReader::import_file(const string file_name) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        import_xyz(file_name);
    else if (file_type == "ckx")
        import_ckx(file_name);
    else if (file_type == "dump")
        import_dump(file_name);
    else
    	require(false, "Unknown file type: " + file_type);

    calc_statistics();
}

const void AtomReader::import_xyz(const string file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type;
    string elem, line, dummy;
    size_t n_atoms;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;
    reserve(n_atoms);

    getline(in_file, line);     // Skip comments line

    // keep storing values from the text file as long as data exists:
    while ((--n_atoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x >> y >> z >> type;
        add_atom(x, y, z, type);
    }
}

const void AtomReader::import_ckx(const string file_name) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    double x, y, z;
    int type;
    string line, dummy;
    size_t n_atoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;

    reserve(n_atoms);

    getline(in_file, line);    // Skip comments line

    // keep storing values from the text file as long as data exists:
    while ((--n_atoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x >> y >> z;
        add_atom(x, y, z, type);
    }
}

const void AtomReader::import_dump(const string file_name) {
    require(false, "AtomReader::import_dump not implemented!");
}

} /* namespace femocs */
