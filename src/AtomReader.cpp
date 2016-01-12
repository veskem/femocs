/*
 * AtomReader.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "AtomReader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cfloat>

using namespace std;
namespace femocs {

AtomReader::AtomReader(const string& fileName) {
    int start = fileName.find_last_of('.') + 1;
	int end = fileName.size();
	string fileType = fileName.substr(start, end);

	if(fileType == "xyz")
		import_xyz(fileName);
	else if(fileType == "ckx")
		import_ckx(fileName);
	else if(fileType == "dump")
		import_dump(fileName);
	else
		cout << "Unknown file type: " << fileType << endl;
}

int AtomReader::import_xyz(const string& fileName) {
	ifstream inFile(fileName, ios::in);
    if (!inFile.is_open()) {
    	cout << "Did not find a file " << fileName << endl;
    	return 1;
    }

    double x,y,z;
    int type;
    string elem, line, dummy;
    size_t nrOfAtoms;
    istringstream iss;

    getline(inFile, line); 	// Read number of atoms
    iss.clear(); iss.str(line);
    iss >> nrOfAtoms;
    this->data.x.reserve(nrOfAtoms);
    this->data.y.reserve(nrOfAtoms);
    this->data.z.reserve(nrOfAtoms);
    this->data.type.reserve(nrOfAtoms);
    this->data.xmin = DBL_MAX; this->data.xmax = DBL_MIN;
    this->data.ymin = DBL_MAX; this->data.ymax = DBL_MIN;
    this->data.zmin = DBL_MAX; this->data.zmax = DBL_MIN;

    getline(inFile, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ( (--nrOfAtoms > 0) && getline(inFile, line) ) {
        iss.clear(); iss.str(line);
		iss >> elem >> x >> y >> z >> type;
    	this->data.x.push_back(x);
        this->data.y.push_back(y);
        this->data.z.push_back(z);
        this->data.type.push_back(type);
        if(x > this->data.xmax) this->data.xmax = x;
        if(y > this->data.ymax) this->data.ymax = y;
        if(z > this->data.zmax) this->data.zmax = z;
        if(x < this->data.xmin) this->data.xmin = x;
        if(y < this->data.ymin) this->data.ymin = y;
        if(z < this->data.zmin) this->data.zmin = z;
    }

    return 0;
}

int AtomReader::import_ckx(const string& fileName) {
	ifstream inFile(fileName, ios::in);
    if (!inFile.is_open()) {
    	cout << "Did not find a file " << fileName << endl;
    	return 1;
    }

    double x,y,z;
    int type;
    string line, dummy;
    size_t nrOfAtoms;
    istringstream iss;

    getline(inFile, line); 	// Read number of atoms
    iss.clear(); iss.str(line);
    iss >> nrOfAtoms;
    this->data.x.reserve(nrOfAtoms);
    this->data.y.reserve(nrOfAtoms);
    this->data.z.reserve(nrOfAtoms);
    this->data.type.reserve(nrOfAtoms);
    this->data.xmin = DBL_MAX; this->data.xmax = DBL_MIN;
    this->data.ymin = DBL_MAX; this->data.ymax = DBL_MIN;
    this->data.zmin = DBL_MAX; this->data.zmax = DBL_MIN;

    getline(inFile, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ( (--nrOfAtoms > 0) && getline(inFile, line) ) {
        iss.clear(); iss.str(line);
		iss >> type >> x >> y >> z;
    	this->data.x.push_back(x);
        this->data.y.push_back(y);
        this->data.z.push_back(z);
        this->data.type.push_back(type);
        if(x > this->data.xmax) this->data.xmax = x;
        if(y > this->data.ymax) this->data.ymax = y;
        if(z > this->data.zmax) this->data.zmax = z;
        if(x < this->data.xmin) this->data.xmin = x;
        if(y < this->data.ymin) this->data.ymin = y;
        if(z < this->data.zmin) this->data.zmin = z;
    }

    return 0;
}

int AtomReader::import_dump(const string& fileName){
	return 1;
}

} /* namespace femocs */
