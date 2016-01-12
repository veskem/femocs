/*
 * AtomReader.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include <vector>
#include <string>

using namespace std;
namespace femocs {

/**
 * Class to import atom coordinates and types
 */
class AtomReader {
public:
	/**
	 * Class to import atomic coordinates and atom types from file
	 * @param file name with type .xyz (Parcas) or .dump (Lammps)
	 */
	AtomReader(const string&);
	virtual ~AtomReader() {};
	struct Data {
	    vector<double> x, y, z;
	    vector<int> type;
	    double xmin, xmax, ymin, ymax, zmin, zmax;
	    const int type_bulk=1, type_surf=2, type_vac=3; 
	};
	Data data;

private:
	/**
	 * Functions to import atoms from file
	 * @param file name
	 * @return error code
	 */
	int import_xyz(const string&);
	int import_ckx(const string&);
	int import_dump(const string&);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
