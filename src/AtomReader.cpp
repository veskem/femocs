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
}

const int AtomReader::get_n_atoms() {
    return this->data.x.size();
}

const void AtomReader::calc_coordination(const Femocs::SimuCell* cell, const double cutoff, const int nnn) {
    if (this->data.simu_type == "md")
        calc_md_coordination(cell, cutoff, nnn);
    else if (this->data.simu_type == "kmc")
        calc_kmc_coordination(cell, nnn);
}

const void AtomReader::calc_md_coordination(const Femocs::SimuCell* cell, const double cutoff, const int nnn) {
    int N = get_n_atoms();
    int i, j, coord;
    double r2, dx, dy, dz;
    double cutoff2 = cutoff*cutoff;

    data.coordination.reserve(N);

//    omp_set_num_threads(this->nrOfOmpThreads);
//#pragma omp parallel for shared(coords,is_surface) private(i,j,dx,dy,dz,r2) reduction(+:coord,N)
    for (i = 0; i < N; ++i) {
        coord = 0;
        for (j = 0; j < N; ++j) {
            if (i == j) continue;
            dx = fabs(data.x[i] - data.x[j]);
            dy = fabs(data.y[i] - data.y[j]);
            dz = fabs(data.z[i] - data.z[j]);

            dx = min(dx, fabs(dx-cell->xbox));  // apply periodic boundary condition in x-direction
            dy = min(dy, fabs(dy-cell->ybox));  // apply periodic boundary condition in y-direction
            //dz = min(dz, fabs(dz-cell->zbox));  // apply periodic boundary condition in z-direction

            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 <= cutoff2) coord++;
            if (coord >= nnn) break; // Coordination can't be bigger than the biggest expected
        }

        data.coordination[i] = coord;
    }
}

const void AtomReader::calc_kmc_coordination(const Femocs::SimuCell* cell, const int nnn) {
    int N = get_n_atoms();
    int i, j, coord;

    data.coordination.reserve(N);

    for (i = 0; i < N; ++i) {
        if(data.type[i] == cell->type_bulk)
            data.coordination[i] = nnn;
        else if(data.type[i] == cell->type_surf)
            data.coordination[i] = (int) nnn/2;
        else
            data.coordination[i] = 0;
    }
}

const void AtomReader::import_helmod() {
    cout << "AtomReader::read_helmod() not implemented!" << endl;
}

const void AtomReader::import_kimocs() {
    cout << "AtomReader::read_kimocs() not implemented!" << endl;
}

const void AtomReader::import_file(const string& file_name, Femocs::SimuCell* cell) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        import_xyz(file_name, cell);
    else if (file_type == "ckx")
        import_ckx(file_name, cell);
    else if (file_type == "dump")
        import_dump(file_name, cell);
    else
        cout << "Unknown file type: " << file_type << endl;

    cell->zmaxbox += cell->zmax;
    cell->xbox = cell->xmax - cell->xmin;
    cell->ybox = cell->ymax - cell->ymin;
    cell->zbox = cell->zmax - cell->zmin;
}

const void AtomReader::import_xyz(const string& file_name, Femocs::SimuCell* cell) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "\nDid not find a file " << file_name << endl;
    }

    double x, y, z;
    int type;
    string elem, line, dummy;
    size_t natoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> natoms;
    init_data(natoms, cell);

    getline(in_file, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ((--natoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x >> y >> z >> type;
        add_data(x, y, z, type, cell);
    }

    return;
}

const void AtomReader::import_ckx(const string& file_name, Femocs::SimuCell* cell) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "Did not find a file " << file_name << endl;
        return;
    }

    double x, y, z;
    int type;
    string line, dummy;
    size_t natoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> natoms;
    init_data(natoms, cell);

    getline(in_file, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ((--natoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x >> y >> z;
        add_data(x, y, z, type, cell);
    }

    return;
}

const void AtomReader::import_dump(const string& file_name, Femocs::SimuCell* cell) {
    cout << "AtomReader::import_dump not implemented!" << endl;
}

const string AtomReader::get_file_type(const string& file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    string file_type = file_name.substr(start, end);

    if (file_type == "xyz") this->data.simu_type = "md";
    if (file_type == "dump") this->data.simu_type = "md";
    if (file_type == "ckx") this->data.simu_type = "kmc";

    return file_type;
}

const void AtomReader::init_data(const int natoms, Femocs::SimuCell* cell) {
    this->data.x.reserve(natoms);
    this->data.y.reserve(natoms);
    this->data.z.reserve(natoms);
    this->data.type.reserve(natoms);
    cell->xmin = DBL_MAX;
    cell->xmax = DBL_MIN;
    cell->ymin = DBL_MAX;
    cell->ymax = DBL_MIN;
    cell->zmin = DBL_MAX;
    cell->zmax = DBL_MIN;
    return;
}

const void AtomReader::add_data(const double x, const double y, const double z, const int type,
        Femocs::SimuCell* cell) {
    this->data.x.push_back(x);
    this->data.y.push_back(y);
    this->data.z.push_back(z);
    this->data.type.push_back(type);

    if (x > cell->xmax) cell->xmax = x;
    if (y > cell->ymax) cell->ymax = y;
    if (z > cell->zmax) cell->zmax = z;
    if (x < cell->xmin) cell->xmin = x;
    if (y < cell->ymin) cell->ymin = y;
    if (z < cell->zmin) cell->zmin = z;
    return;
}

} /* namespace femocs */
