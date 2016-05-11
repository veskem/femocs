/*
 * Medium.h
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

#include "Macros.h"
#include "Vec3.h"
#include <deal.II/numerics/vector_tools.h>

using namespace std;
namespace femocs {

/** Class to define elementary operations between points with double coordinates*/
class Point3d{
public:

    Point3d(double xx, double yy, double zz) :
            x(xx), y(yy), z(zz) {}

    double distance(const Point3d &n) {
        Vec3d dif(x - n.x, y - n.y, z - n.z);
        return dif.length();
    }

    double distance(const dealii::Point<3> &p) {
        Vec3d dif(x - p[0], y - p[1], z - p[2]);
        return dif.length();
    }

    bool operator ==(const dealii::Point<3> &p) {
        return x == p[0] && y == p[1] && z == p[2];
    }

    // Point3d accessors
    const double& operator [](uint8_t i) const {
        return (&x)[i];
    }
    double& operator [](uint8_t i) {
        return (&x)[i];
    }

    double x,y,z;
};

class Medium {
public:

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
    const void output(const string file_name);

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

    /** Return x-coordinate of i-th atom */
    const double get_x(const int i);
    /** Return y-coordinate of i-th atom */
    const double get_y(const int i);
    /** Return z-coordinate of i-th atom */
    const double get_z(const int i);
    /** Return coordination of i-th atom */
    const int get_coordination(const int i);

    /** Return number of atoms in a Medium */
    const int get_n_atoms();

    struct Sizes {
        double xmin;    //!< Minimum value of x-coordinate
        double xmax;    //!< Maximum value of x-coordinate
        double ymin;    //!< Minimum value of y-coordinate
        double ymax;    //!< Maximum value of y-coordinate
        double zmin;    //!< Minimum value of z-coordinate
        double zmax;    //!< Maximum value of z-coordinate
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

    /** Get i-th entry from all data vectors */
    const string get_data_string(const int i);
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
