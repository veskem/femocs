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
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;
namespace femocs {

AtomReader::AtomReader() {
    init_statistics();
}

const void AtomReader::init_statistics() {
    sizes.xmin = sizes.ymin = sizes.zmin = DBL_MAX;
    sizes.xmax = sizes.ymax = sizes.zmax = DBL_MIN;
}

// Calculate the statistics about coordinates in AtomReader
const void AtomReader::calc_statistics() {
    init_statistics();

    for (int i = 0; i < get_n_atoms(); ++i) {
        if (sizes.xmax < get_x(i)) sizes.xmax = get_x(i);
        if (sizes.xmin > get_x(i)) sizes.xmin = get_x(i);
        if (sizes.ymax < get_y(i)) sizes.ymax = get_y(i);
        if (sizes.ymin > get_y(i)) sizes.ymin = get_y(i);
        if (sizes.zmax < get_z(i)) sizes.zmax = get_z(i);
        if (sizes.zmin > get_z(i)) sizes.zmin = get_z(i);
    }

    sizes.xbox = sizes.xmax - sizes.xmin;
    sizes.ybox = sizes.ymax - sizes.ymin;
    sizes.zbox = sizes.zmax - sizes.zmin;
    sizes.zminbox = sizes.zmin;
    sizes.zmaxbox = sizes.zmax;
}

const void AtomReader::resize_box(const double zmin, const double zmax) {
    sizes.zminbox = zmin;
    sizes.zmaxbox = zmax;
    sizes.zbox = zmax - zmin;
}

const int AtomReader::get_n_atoms() {
    return x.size();
}

const double AtomReader::get_x(const int i) {
    return x[i];
}

const double AtomReader::get_y(const int i) {
    return y[i];
}

const double AtomReader::get_z(const int i) {
    return z[i];
}

const int AtomReader::get_type(const int i) {
    return type[i];
}

const int AtomReader::get_coord(const int i) {
    return coordination[i];
}

const void AtomReader::calc_coordination(const double cutoff, const int nnn) {
    if (types.simu_type == "md")
        calc_md_coordination(cutoff, nnn);
    else if (types.simu_type == "kmc") calc_kmc_coordination(nnn);
}

const void AtomReader::calc_md_coordination(const double cutoff, const int nnn) {
    int N = get_n_atoms();
    int i, j, coord;
    double r2, dx, dy, dz;
    double cutoff2 = cutoff * cutoff;

    coordination.reserve(N);

//    omp_set_num_threads(this->nrOfOmpThreads);
//#pragma omp parallel for shared(coords,is_surface) private(i,j,dx,dy,dz,r2) reduction(+:coord,N)
    for (i = 0; i < N; ++i) {
        coord = 0;
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

const void AtomReader::calc_kmc_coordination(const int nnn) {
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

const void AtomReader::import_helmod() {
    cout << "AtomReader::read_helmod() not implemented!" << endl;
}

const void AtomReader::import_kimocs() {
    cout << "AtomReader::read_kimocs() not implemented!" << endl;
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
        cout << "Unknown file type: " << file_type << endl;

    calc_statistics();
}

const void AtomReader::import_xyz(const string file_name) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "\nDid not find a file " << file_name << endl;
    }

    double x, y, z;
    int type;
    string elem, line, dummy;
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
        iss >> elem >> x >> y >> z >> type;
        add_atom(x, y, z, type);
    }
}

const void AtomReader::import_ckx(const string file_name) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "\nDid not find a file " << file_name << endl;
    }

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
    cout << "AtomReader::import_dump not implemented!" << endl;
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

const void AtomReader::reserve(const int n_atoms) {
    x.reserve(n_atoms);
    y.reserve(n_atoms);
    z.reserve(n_atoms);
    type.reserve(n_atoms);
}

const void AtomReader::add_atom(const double x, const double y, const double z, const int type) {
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->type.push_back(type);
}

} /* namespace femocs */
