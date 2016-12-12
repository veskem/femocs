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

/** Class to import atoms from atomistic simulation and to divide them into different categories */
class AtomReader: public Medium {
public:
    /** Constructor for AtomReader */
    AtomReader();

    /** Get type of i-th atom in AtomReader */
    const int get_type(const int i) const;

    /**
     * Function to import file with atom coordinates and types
     * @param file_name  path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     */
    const void import_file(const string file_name);

    /** Function to transform atomic data from Parcas format into AtomReader one */
    const void import_parcas(int n_atoms, const double* coordinates, const double* box);

    /** Function to transform atomic data from Helmod format into AtomReader one */
    const void import_helmod(int n_atoms, double* x, double* y, double* z, int* types);

    /** Function to transform atomic data from Kimocs format into AtomReader one */
    const void import_kimocs();

    /**
     * Calculate coordination for all the atoms in AtomReader
     * @param cutoff   cut off radius for coordination analysis
     * @param nnn      number of nearest neighbours in a crystal
     * @param nborlist neighbour list of atoms
     */
    const void calc_coordination(const int nnn, const double cutoff, const int* nborlist);

    /** Calculate coordination for all the atoms using brute force technique */
    const void calc_coordination(const double cutoff);

    /** Calculate pseudo-coordination for all the atoms using the atom types */
    const void calc_coordination(const int nnn);

    /** Extract atom types from calculated atom coordinations */
    const void extract_types(const int nnn, const double latconst);

    /** Redefine the min and max values for z-coordinates */
    const void resize_box(const double zmin, const double zmax);

    /** Compare current and previous run and detect whether 
     * any of the atoms has moved more than threshold */
    const bool equals_previous_run(const double eps);
    
    /** Calculate the root mean average distance the atoms have moved 
     * between previous and current run */
    const double diff_from_prev_run(const double eps);
    
    /** Store the atom coordinates from current run */
    const void save_current_run_points(const double eps);

private:
    vector<int> types;               ///< types of atoms
    vector<Point3> previous_point;   ///< atom coordinates from previous run

    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    const void import_xyz(const string file_name);
    const void import_ckx(const string file_name);
    const void import_dump(const string file_name);

    /** Function to detect evaporated atoms by their coordinations;
     * an atom is considered evaporated if it's coordiantion is between the defined limits. */
    const void check_coordination();

    /** Reserve memory for data vectors */
    const void reserve(const int n_atoms);
    
    /** Add atom with its id and type */
    const void add_atom(const int id, const Point3 &point, const int type);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
