/*
 * Medium.h
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

#include "Macros.h"
#include "Primitives.h"

using namespace std;
namespace femocs {

class Medium {
public:
    /** Medium constructor */
    Medium();
    virtual ~Medium() {}; // = 0;

    /** Sort the atoms by their radial coordinate from origin */
    const void sort_atoms(const Point2 &origin, const string& direction = "up");

    /** Define the addition of two Mediums */
    Medium& operator +=(Medium &m);

    /** Add data from other Medium to current one */
    const void add(Medium *m);

    /** Reserve memory for data vectors */
    virtual const void reserve(const int n_atoms);

    /** Add Atom to the system */
    const void add_atom(const Atom& atom);

    /**
     * Export the data of Medium to file
     * @param file_name - path for file to save the data
     */
    void output(const string file_name);

    /** Calculate statistics about the coordinates in Medium */
    const void calc_statistics();
    /** Set the id of i-th atom */
    const void set_id(const int i, const int id);
    /** Set the coordinates of i-th atom */
    const void set_point(const int i, const Point3& p);
    /** Set the x-coordinate of i-th atom */
    const void set_x(const int i, const double x);
    /** Set the y-coordinate of i-th atom */
    const void set_y(const int i, const double y);
    /** Set the z-coordinate of i-th atom */
    const void set_z(const int i, const double z);
    /** Set the coordination of i-th atom */
    const void set_coordination(const int i, const int coord);

    /** Return i-th Atom */
    const Atom get_atom(const int i);
    /** Return x-, y- and z-coordinate and associated operators for i-th atom */
    const Point3 get_point(const int i);
    /** Return x- and y-coordinate and associated operators for i-th atom */
    const Point2 get_point2(const int i);
    /** Return ID of i-th atom */
    const int get_id(const int i);
    /** Return coordination of i-th atom */
    const int get_coordination(const int i);
    /** Return number of atoms in a Medium */
    const int get_n_atoms();

    /** Statistics about system size */
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
        double xmean;   //!< average value of x-coordinate
        double ymean;   //!< average value of y-coordinate
        double zmean;   //!< average value of z-coordinate
    } sizes;

    /** Statistics about crystal structure */
    struct CrysStruct {
        double latconst;    //!< Lattice constant
        int nnn;            //!< Number of nearest neighbours
    } crys_struct;

protected:
    vector<Atom> atoms;

    /** Initialise statistics about the coordinates in Medium */
    const void init_statistics();

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    virtual const string get_data_string(const int i);
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
