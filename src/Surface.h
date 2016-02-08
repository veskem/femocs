/*
 * Medium.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <string>
#include <vector>

using namespace std;

namespace femocs {

/**
 * Routines and data related to making a surface
 */
class Surface {
public:
    /**
     * Surface constructor
     * @param nr_of_atoms - number of atoms on a surface
     */
    Surface(const int natoms, const double latconst);
    ~Surface() {
    }
    ;

    /**
     * Function to add atom to Surface
     * @param x - x coordinate of atom
     * @param y - y coordinate of atom
     * @param z - z coordinate of atom
     * @param coord - coordination of atom
     * @param type - atom type
     */
    const void add(const double x, const double y, const double z, const int coord, const int type);

    /**
     * Function to export the data of Surface.
     * @param file_name - path for file to save the data
     */
    const void output(const string& file_name);

    vector<double>* getXs();
    vector<double>* getYs();
    vector<double>* getZs();

    /** Return x-coordinate of i-th atom */
    const double getX(const int i);
    /** Return y-coordinate of i-th atom */
    const double getY(const int i);
    /** Return z-coordinate of i-th atom */
    const double getZ(const int i);
    /** Return type of i-th atom */
    const int getType(const int i);
    /** Return number of atoms in a Surface */
    const int getN();
    /** Return lattice constant of extracted material */
    const double getLatconst();

private:
    vector<double> x, y, z;	//!< Real coordinates
    vector<double> Sx, Sy, Sz; //!< Surface area components
    vector<int> coordination;	//!< Number of nearest neighbours in cutoff radius
    vector<int> type;			//!< Atom type (-1-fixed, 1-bulk, 2-surface, 3-vacancy)
    vector<bool> isEvaporated; //!< List of evaporated atoms
    double latconst;            //!< lattice constant of extracted material
};

} /* namespace femocs */

#endif /* SURFACE_H_ */
