/*
 * Medium.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Surface.h"

#include <fstream>

using namespace std;

namespace femocs {

// Constructor initializes arrays
Surface::Surface(const int natoms, double latconst) {
    this->x.reserve(natoms);				// Default values are 0
    this->y.reserve(natoms);
    this->z.reserve(natoms);
    this->Sx.reserve(natoms);
    this->Sy.reserve(natoms);
    this->Sz.reserve(natoms);
    this->coordination.reserve(natoms);
    this->type.reserve(natoms);
    this->isEvaporated.reserve(natoms);	// Default values are False
    this->latconst = latconst;
}
;

// Add atom coordinates
const void Surface::add(const double x, const double y, const double z, const int coord,
        const int type) {
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->coordination.push_back(coord);
    this->type.push_back(type);
}

// Get element from x, y or z vector
const double Surface::getX(const int i) {
    return x[i];
}
const double Surface::getY(const int i) {
    return y[i];
}
const double Surface::getZ(const int i) {
    return z[i];
}
// Get lattice constant
const double Surface::getLatconst() {
    return latconst;
}
// Get atom type
const int Surface::getType(const int i) {
    return type[i];
}

// Get number of atoms on surface
const int Surface::getN() {
    return x.size();
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
