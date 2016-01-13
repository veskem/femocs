/*
 * Medium.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Medium.h"

#include <fstream>

using namespace std;

namespace femocs {

// Constructor initializes arrays
Surface::Surface(int nr_of_atoms) {
    this->x.reserve(nr_of_atoms);				// Default values are 0
    this->y.reserve(nr_of_atoms);
    this->z.reserve(nr_of_atoms);
    this->Sx.reserve(nr_of_atoms);
    this->Sy.reserve(nr_of_atoms);
    this->Sz.reserve(nr_of_atoms);
    this->coordination.reserve(nr_of_atoms);
    this->type.reserve(nr_of_atoms);
    this->isEvaporated.reserve(nr_of_atoms);	// Default values are False
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
