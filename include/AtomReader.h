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
    int get_type(const int i) const;

    /**
     * Function to import file with atom coordinates and types
     * @param file_name  path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     */
    void import_file(const string &file_name);

    /** Function to transform atomic data from Parcas format into AtomReader one */
    void import_parcas(const int n_atoms, const double* coordinates, const double* box);

    /** Function to transform atomic data from Helmod format into AtomReader one */
    void import_helmod(const int n_atoms, const double* x, const double* y, const double* z, const int* types);

    /** Function to transform atomic data from Kimocs format into AtomReader one */
    void import_kimocs();

    /**
     * Calculate coordination for all the atoms in AtomReader
     * @param cutoff   cut off radius for coordination analysis
     * @param nnn      number of nearest neighbours in a crystal
     * @param nborlist neighbour list of atoms
     */
    void calc_coordination(const int nnn, const double cutoff, const int* nborlist);

    /** Calculate coordination for all the atoms using brute force technique */
    void calc_coordination(const double cutoff);

    /** Calculate pseudo-coordination for all the atoms using the atom types */
    void calc_coordination(const int nnn);

    /** Function to detect evaporated atoms by their coordinations;
     * an atom is considered evaporated if it's coordiantion is between the defined limits. */
    void check_coordination();

    /** Extract atom types from calculated atom coordinations */
    void extract_types(const int nnn, const double latconst);

    /** Redefine the min and max values for z-coordinates */
    void resize_box(const double zmin, const double zmax);
    
    /** Redefine the min and max values for x, y and z - coordinates */
    void resize_box(const double xmin, const double xmax, const double ymin, const double ymax,
            const double zmin, const double zmax);

    /** Compare current and previous run and detect whether 
     * any of the atoms has moved more than threshold */
    bool equals_previous_run(const double eps);
    
    /** Calculate the root mean average distance the atoms have moved 
     * between previous and current run */
    double get_rms_distance(const double eps);
    
    /** Store the atom coordinates from current run */
    void save_current_run_points(const double eps);

    double rms_distance;           ///< rms distance between atoms from previous and current run
private:
    vector<int> types;             ///< types of atoms
    vector<Point3> previous_point; ///< atom coordinates from previous run

    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    void import_xyz(const string &file_name);
    void import_ckx(const string &file_name);
    void import_dump(const string &file_name);

    /** Reserve memory for data vectors */
    void reserve(const int n_atoms);
    
    /** Add atom with its id and type */
    void add_atom(const int id, const Point3 &point, const int type);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
