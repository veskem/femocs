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

#include "Femocs.h"

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
    const void import_file(const string& file_name, Femocs::SimuCell* cell);

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

    /** Get coordinatino of i-th atom in AtomReader */
    const int get_coord(const int i);

    /**
     * Calculate coordination for all the atoms in AtomReader
     * @param cell - data where simulation cell sizes is saved
     * @param cutoff - cut off radius for coordination analysis
     * @param nnn - number of nearest neighbours in a crystal
     */
    const void calc_coordination(const Femocs::SimuCell* cell, const double cutoff, const int nnn);

    /** Struct for holding data of the whole simulation cell. */
    struct Data {
        double xmin;        //!< minimum x-coordinate in simulation cell
        double xmax;        //!< maximum x-coordinate in simulation cell
        double ymin;        //!< minimum y-coordinate in simulation cell
        double ymax;        //!< maximum y-coordinate in simulation cell
        double zmin;        //!< minimum z-coordinate in simulation cell
        double zmax;        //!< maximum z-coordinate in simulation cell
        string simu_type;   //!< Type of simulation; md or kmc
    };
    Data data; //!< AtomReader data

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
    const string get_file_type(const string& file_name);
    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    const void import_xyz(const string& file_name, Femocs::SimuCell* cell);
    const void import_ckx(const string& file_name, Femocs::SimuCell* cell);
    const void import_dump(const string& file_name, Femocs::SimuCell* cell);

    const void calc_md_coordination(const Femocs::SimuCell* cell, const double cutoff, const int nnn);
    const void calc_kmc_coordination(const Femocs::SimuCell* cell, const int nnn);

    const void init_data(const int natoms, Femocs::SimuCell* cell);
    const void add_atom(const double x, const double y, const double z, const int type,
            Femocs::SimuCell* cell);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
