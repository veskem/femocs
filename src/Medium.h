/*
 * Medium.h
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <float.h>

using namespace std;
namespace femocs {

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
    const void output(const string& file_name);
    
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
