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
    const void add_atom(const double x, const double y, const double z, const int coord, const int type);

    /**
     * Function to export the data of Surface.
     * @param file_name - path for file to save the data
     */
    const void output(const string& file_name);

    vector<double>* get_xs();
    vector<double>* get_ys();
    vector<double>* get_zs();

    /** Set the x-coordinate of i-th atom */
    void set_x(const int i, const double x);
    /** Set the y-coordinate of i-th atom */
    void set_y(const int i, const double y);
    /** Set the z-coordinate of i-th atom */
    void set_z(const int i, const double z);
    /** Set the type of i-th atom */
    void set_type(const int i, const int type);
    /** Set the coordination of i-th atom */
    void set_coordination(const int i, const int coord);

    /** Reserve N entries to the Surface data arrays */
    const void reserve(const int N);

    /** Return x-coordinate of i-th atom */
    const double get_x(const int i);
    /** Return y-coordinate of i-th atom */
    const double get_y(const int i);
    /** Return z-coordinate of i-th atom */
    const double get_z(const int i);
    /** Return type of i-th atom */
    const int get_type(const int i);
    /** Return coordination of i-th atom */
    const int get_coordination(const int i);
    /** Return number of atoms in a Surface */
    const int get_n_atoms();
    /** Calculate statistics about Surface atoms */
    const void calc_statistics();

    struct Sizes {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        double latconst;
    };

    Sizes sizes;

private:
    vector<double> x, y, z;	//!< Real coordinates
    vector<double> Sx, Sy, Sz; //!< Surface area components
    vector<int> coordination;	//!< Number of nearest neighbours in cutoff radius
    vector<int> type;			//!< Atom type (-1-fixed, 1-bulk, 2-surface, 3-vacancy)
    vector<bool> isEvaporated; //!< List of evaporated atoms
//    double latconst;            //!< lattice constant of extracted material

    const void init_statistics();
};

} /* namespace femocs */

#endif /* SURFACE_H_ */
