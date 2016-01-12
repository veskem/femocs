/*
 * Medium.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Medium.h"
#include "AtomReader.h"
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

namespace femocs {

// Constructor initializes arrays
Surface::Surface(int nrOfAtoms) {
	this->x.reserve(nrOfAtoms);				// Default values are 0
	this->y.reserve(nrOfAtoms);
	this->z.reserve(nrOfAtoms);
	this->Sx.reserve(nrOfAtoms);
	this->Sy.reserve(nrOfAtoms);
	this->Sz.reserve(nrOfAtoms);
	this->coordination.reserve(nrOfAtoms);
	this->type.reserve(nrOfAtoms);
	this->isEvaporated.reserve(nrOfAtoms);	// Default values are False
};

// Add atom coordinates
const void Surface::add(double x, double y, double z, int coord, int type) {
	this->x.push_back(x);
	this->y.push_back(y);
	this->z.push_back(z);
	this->coordination.push_back(coord);
	this->type.push_back(type);
}

// Get element from x, y or z vector
const double Surface::getX(const int i) { return x[i]; }
const double Surface::getY(const int i) { return y[i]; }
const double Surface::getZ(const int i) { return z[i]; }

// Get atom type
const int Surface::getType(const int i) { return type[i]; }

// Get number of atoms on surface
const int Surface::getN(){ return x.size(); }


// Output surface data to file
const void Surface::output(const string& fileName) {
	int nrOfAtoms = this->x.size();
	ofstream myfile;
	myfile.open (fileName);
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

//void Surface::smooth(Smoother s) {}

} /* namespace femocs */
