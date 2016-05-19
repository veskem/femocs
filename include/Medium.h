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

    /** Define the addition of two Mediums */
    Medium& operator +=(Medium &m) {
        int n_atoms1 = get_n_atoms();
        int n_atoms2 = m.get_n_atoms();

        this->reserve(n_atoms1 + n_atoms2);

        for(int i = 0; i < n_atoms2; ++i)
            add_atom(m.get_x(i), m.get_y(i), m.get_z(i), m.get_coordination(i));

        this->calc_statistics();
        return *this;
    }

    const void add(Medium *m) {
        int n_atoms1 = get_n_atoms();
        int n_atoms2 = m->get_n_atoms();

        this->reserve(n_atoms1 + n_atoms2);

        for(int i = 0; i < n_atoms2; ++i)
            add_atom(m->get_x(i), m->get_y(i), m->get_z(i), m->get_coordination(i));

        this->calc_statistics();
    }

    /** Reserve memory for data vectors */
    const void reserve(const int n_atoms);

    /**
     * Add atom with its parameters to the data vectors
     * @param x - x-coordinate of the atom
     * @param y - y-coordinate of the atom
     * @param z - z-coordinate of the atom
     * @param coord - coordination of the atom; 0 in case of none
     */
    const void add_atom(const double x, const double y, const double z, const int coord);

    /**
     * Function to export the data of Medium
     * @param file_name - path for file to save the data
     */
//    const void output(const string file_name);

    /** Calculate statistics about the coordinates in Medium */
    const void calc_statistics();

    /** Set the x-coordinate of i-th atom */
    void set_x(const int i, const double x);
    /** Set the y-coordinate of i-th atom */
    void set_y(const int i, const double y);
    /** Set the z-coordinate of i-th atom */
    void set_z(const int i, const double z);
    /** Set the coordination of i-th atom */
    void set_coordination(const int i, const int coord);

    /** Return x-, y- and z-coordinate and associated operators for i-th atom */
    Point3d get_point(const int i);

    /** Return x- and y-coordinate and associated operators for i-th atom */
    const Point2d get_point2d(const int i);

    /** Return x-coordinate of i-th atom */
    const double get_x(const int i);
    /** Return y-coordinate of i-th atom */
    const double get_y(const int i);
    /** Return z-coordinate of i-th atom */
    const double get_z(const int i);
    /** Return coordination of i-th atom */
    const int get_coordination(const int i);

    /** Return number of atoms in a Medium */
    int get_n_atoms();

    /** Function to export the data to file */
    const void output(const string file_name);

    struct Sizes {
        double xmin;    //!< Minimum value of x-coordinate
        double xmax;    //!< Maximum value of x-coordinate
        double xmean;   //!< Average value of x-coordinate
        double ymin;    //!< Minimum value of y-coordinate
        double ymax;    //!< Maximum value of y-coordinate
        double ymean;   //!< Average value of y-coordinate
        double zmin;    //!< Minimum value of z-coordinate
        double zmax;    //!< Maximum value of z-coordinate
        double zmean;   //!< Average value of z-coordinate
    };

    /** Statistics about system size */
    Sizes sizes;

    struct CrysStruct {
        double latconst;    //!< Lattice constant
        int nnn;            //!< Number of nearest neighbours
    };

    /** Statistics about crystal structure */
    CrysStruct crys_struct;

protected:
    vector<double> x, y, z;   //!< Atom coordinates
    vector<int> coordination; //!< Atom coordination - nr of nearest neighbours within cut off radius

    /** Initialise statistics about the coordinates in Medium */
    const void init_statistics();

    /** Function to extract file extension from file name */
    const string get_file_type(const string file_name);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
