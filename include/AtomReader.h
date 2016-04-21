/*
 * AtomReader.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include "Macros.h"

using namespace std;
namespace femocs {

/**
 * Class to import atom coordinates and types.
 */
class AtomReader {
public:
    /**
     * Constructor for AtomReader.
     */
    AtomReader();
    virtual ~AtomReader() {
    }
    ;

    /**
     * Function to import file with atom coordinates and types.
     * @param file_name - path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     * @param cell - pointer to struct holding the data about simulation cell size
     */
    const void import_file(const string file_name);

    /**
     * Function to transform atomic data from Helmod format into AtomReader one
     */
    const void import_helmod();
    /**
     * Function to transform atomic data from Kimocs format into AtomReader one
     */
    const void import_kimocs();

    /** Get number of atoms in AtomReader */
    const int get_n_atoms();

    /** Get x-coordinate of i-th atom in AtomReader */
    const double get_x(const int i);

    /** Get y-coordinate of i-th atom in AtomReader */
    const double get_y(const int i);

    /** Get z-coordinate of i-th atom in AtomReader */
    const double get_z(const int i);

    /** Get type of i-th atom in AtomReader */
    const int get_type(const int i);

    /** Get coordination of i-th atom in AtomReader */
    const int get_coord(const int i);

    /**
     * Calculate coordination for all the atoms in AtomReader
     * @param cell - data where simulation cell sizes is saved
     * @param cutoff - cut off radius for coordination analysis
     * @param nnn - number of nearest neighbours in a crystal
     */
    const void calc_coordination(const double cutoff, const int nnn);

    /** Calculate statistics about the coordinates in AtomReader */
    const void calc_statistics();

    /** Redefine the min and max values for z-coordinates */
    const void resize_box(const double zmin, const double zmax);

    /** Struct for holding data of the whole simulation cell. */
    struct Types {
        int type_bulk = 1;  //!< type of bulk material
        int type_surf = 2;  //!< type of open material surface
        int type_vacancy = 3; //!< type of vacancies
        int type_vacuum = 3;  //!< type of vacuum
        int type_edge = 0;  //!< type of the rim/outer edge of surface
        int type_fixed = -1;  //!< type of fixed atoms
        int type_xmin = 4;  //!< type of atom on negative x-face of simulation cell
        int type_ymin = 5;  //!< type of atom on negative y-face of simulation cell
        int type_zmin = 6;  //!< type of atom on negative z-face of simulation cell
        int type_xmax = 10; //!< type of atom on positive x-face of simulation cell
        int type_ymax = 9;  //!< type of atom on positive y-face of simulation cell
        int type_zmax = 8;  //!< type of atom on positive z-face of simulation cell
        int type_none = 7;  //!< type of atom with unknown position

        string simu_type;   //!< Type of simulation; md | kmc | standalone
    };

    struct Sizes {
        double xmin;    //!< minimum x-coordinate of atoms
        double xmax;    //!< maximum x-coordinate of atoms
        double ymin;    //!< minimum y-coordinate of atoms
        double ymax;    //!< maximum y-coordinate of atoms
        double zmin;    //!< minimum z-coordinate of atoms
        double zmax;    //!< maximum z-coordinate of atoms
        double zminbox; //!< minimum z-coordinate of simulation box
        double zmaxbox; //!< maximum z-coordinate of simulation box
        double xbox;    //!< simulation box size in x-direction
        double ybox;    //!< simulation box size in y-direction
        double zbox;    //!< simulation box size in z-direction
    };

    /** Types of regions used in the simulation */
    Types types;

    /** Statistics about system size */
    Sizes sizes;

private:
    vector<double> x;   //!< x-coordinates of atoms
    vector<double> y;   //!< y-coordinates of atoms
    vector<double> z;   //!< z-coordinates of atoms
    vector<int> type;   //!< types of atoms
    vector<int> coordination;  //!< coordinations of atoms

    /**
     * Function to extract file extension from file name
     * @param file_name - name of the file with extension
     */
    const string get_file_type(const string file_name);

    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    const void import_xyz(const string file_name);
    const void import_ckx(const string file_name);
    const void import_dump(const string file_name);

    const void calc_md_coordination(const double cutoff, const int nnn);
    const void calc_kmc_coordination(const int nnn);

    const void reserve(const int n_atoms);
    const void add_atom(const double x, const double y, const double z, const int type);

    /** Initialise statistics about coordinates in AtomReader */
    const void init_statistics();
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */