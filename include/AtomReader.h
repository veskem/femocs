/*
 * AtomReader.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include "Macros.h"
#include "Medium.h"

using namespace std;
namespace femocs {

/** Class to import atom coordinates and types. */
class AtomReader: public Medium {
public:
    /** Constructor for AtomReader. */
    AtomReader();

    /** Get type of i-th atom in AtomReader */
    const int get_type(const int i);

    /**
     * Function to import file with atom coordinates and types.
     * @param file_name - path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     * @param cell - pointer to struct holding the data about simulation cell size
     */
    const void import_file(const string file_name);

    /** Function to transform atomic data from Parcas format into AtomReader one */
    const void import_parcas(int n_atoms, const double* coordinates, const double* box);

    /** Function to transform atomic data from Helmod format into AtomReader one */
    const void import_helmod(int n_atoms, double* x, double* y, double* z, int* types);

    /** Function to transform atomic data from Kimocs format into AtomReader one */
    const void import_kimocs();

    const void output(const string file_name);

    /**
     * Calculate coordination for all the atoms in AtomReader
     * @param cell - data where simulation cell sizes is saved
     * @param cutoff - cut off radius for coordination analysis
     * @param nnn - number of nearest neighbours in a crystal
     */
    const void calc_coordination(int nnn, double cutoff = 0.0);

    const void calc_coordination(int nnn, double cutoff, const int* nborlist);

    const void extract_types(int nnn);

//    /** Calculate statistics about the coordinates in AtomReader */
//    const void calc_statistics();

    /** Redefine the min and max values for z-coordinates */
    const void resize_box(const double zmin, const double zmax);

//    /** Types of regions used in the simulation */
//    struct Types {
//        int TYPE_BULK = 1;    //!< type of bulk material
//        int TYPE_SURFACE = 2; //!< type of open material surface
//        int TYPE_VACANCY = 3; //!< type of vacancies
//        int TYPE_VACUUM = 3;  //!< type of vacuum
//        int TYPE_EDGE = 0;    //!< type of the rim/outer edge of surface
//        int TYPE_FIXED = -1;  //!< type of fixed atoms
//        int TYPE_XMIN = 4;  //!< type of atom on negative x-face of simulation cell
//        int TYPE_YMIN = 5;  //!< type of atom on negative y-face of simulation cell
//        int TYPE_ZMIN = 6;  //!< type of atom on negative z-face of simulation cell
//        int TYPE_XMAX = 10; //!< type of atom on positive x-face of simulation cell
//        int TYPE_YMAX = 9;  //!< type of atom on positive y-face of simulation cell
//        int TYPE_ZMAX = 8;  //!< type of atom on positive z-face of simulation cell
//        int TYPE_NONE = 7;  //!< type of atom with unknown position
//    } types;

private:
    vector<int> type;   //!< types of atoms

    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    const void import_xyz(const string file_name);
    const void import_ckx(const string file_name);
    const void import_dump(const string file_name);

    const void calc_slow_coordination(const double cutoff, const int nnn);
    const void calc_dummy_coordination(const int nnn);

    const void reserve(const int n_atoms);
    const void add_atom(const int id, const Point3 &point, const int type);

//    /** Initialise statistics about coordinates in AtomReader */
//    const void init_statistics();

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
