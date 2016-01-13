/*
 * AtomReader.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include <string>
#include <vector>

using namespace std;
namespace femocs {

/**
 * Class to import atom coordinates and types.
 */
class AtomReader {
public:
    /**
     * Constructor for AtomReader.
     * @param file_name - path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     */
    AtomReader(const string& file_name);
    virtual ~AtomReader() {};

    /** Struct for holding data of the whole simulation cell. */
    struct Data {
        vector<double> x;   //!< x-coordinates of atoms
        vector<double> y;   //!< y-coordinates of atoms
        vector<double> z;   //!< z-coordinates of atoms
        vector<int> type;   //!< types of atoms
        double xmin;        //!< minimum x-coordinate in simulation cell
        double xmax;        //!< maximum x-coordinate in simulation cell
        double ymin;        //!< minimum y-coordinate in simulation cell
        double ymax;        //!< maximum y-coordinate in simulation cell
        double zmin;        //!< minimum z-coordinate in simulation cell
        double zmax;        //!< maximum z-coordinate in simulation cell
//	    const int type_bulk=1, type_surf=2, type_vac=3;
    };
    Data data; //!< AtomReader data

private:
    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    const void import_xyz(const string& file_name);
    const void import_ckx(const string& file_name);
    const void import_dump(const string& file_name);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
