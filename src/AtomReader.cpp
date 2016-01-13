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

using namespace std;
namespace femocs {

AtomReader::AtomReader(const string& file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    string fileType = file_name.substr(start, end);

    if (fileType == "xyz")
        import_xyz(file_name);
    else if (fileType == "ckx")
        import_ckx(file_name);
    else if (fileType == "dump")
        import_dump(file_name);
    else
        cout << "Unknown file type: " << fileType << endl;
}

const void AtomReader::import_xyz(const string& file_name) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "Did not find a file " << file_name << endl;
        return;
    }

    double x, y, z;
    int type;
    string elem, line, dummy;
    size_t nrOfAtoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> nrOfAtoms;
    this->data.x.reserve(nrOfAtoms);
    this->data.y.reserve(nrOfAtoms);
    this->data.z.reserve(nrOfAtoms);
    this->data.type.reserve(nrOfAtoms);
    this->data.xmin = DBL_MAX;
    this->data.xmax = DBL_MIN;
    this->data.ymin = DBL_MAX;
    this->data.ymax = DBL_MIN;
    this->data.zmin = DBL_MAX;
    this->data.zmax = DBL_MIN;

    getline(in_file, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ((--nrOfAtoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x >> y >> z >> type;
        this->data.x.push_back(x);
        this->data.y.push_back(y);
        this->data.z.push_back(z);
        this->data.type.push_back(type);
        if (x > this->data.xmax) this->data.xmax = x;
        if (y > this->data.ymax) this->data.ymax = y;
        if (z > this->data.zmax) this->data.zmax = z;
        if (x < this->data.xmin) this->data.xmin = x;
        if (y < this->data.ymin) this->data.ymin = y;
        if (z < this->data.zmin) this->data.zmin = z;
    }

    return;
}

const void AtomReader::import_ckx(const string& file_name) {
    ifstream in_file(file_name, ios::in);
    if (!in_file.is_open()) {
        cout << "Did not find a file " << file_name << endl;
        return;
    }

    double x, y, z;
    int type;
    string line, dummy;
    size_t nrOfAtoms;
    istringstream iss;

    getline(in_file, line); 	// Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> nrOfAtoms;
    this->data.x.reserve(nrOfAtoms);
    this->data.y.reserve(nrOfAtoms);
    this->data.z.reserve(nrOfAtoms);
    this->data.type.reserve(nrOfAtoms);
    this->data.xmin = DBL_MAX;
    this->data.xmax = DBL_MIN;
    this->data.ymin = DBL_MAX;
    this->data.ymax = DBL_MIN;
    this->data.zmin = DBL_MAX;
    this->data.zmax = DBL_MIN;

    getline(in_file, line);    // Skip comments line

    // keep storing values from the text file so long as data exists:
    while ((--nrOfAtoms > 0) && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x >> y >> z;
        this->data.x.push_back(x);
        this->data.y.push_back(y);
        this->data.z.push_back(z);
        this->data.type.push_back(type);
        if (x > this->data.xmax) this->data.xmax = x;
        if (y > this->data.ymax) this->data.ymax = y;
        if (z > this->data.zmax) this->data.zmax = z;
        if (x < this->data.xmin) this->data.xmin = x;
        if (y < this->data.ymin) this->data.ymin = y;
        if (z < this->data.zmin) this->data.zmin = z;
    }

    return;
}

const void AtomReader::import_dump(const string& file_name) {
    return;
}

} /* namespace femocs */
