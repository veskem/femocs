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
    Medium(const int n_atoms);
    virtual ~Medium() {}; // = 0;

    /** Sort the atoms by their x, y or z coordinate (coord=0|1|2) or radial coordinate(coord=3) */
    const void sort_atoms(const int coord, const string& direction="up");

    /** Sort the atoms twice by their x, y or z coordinate */
    const void sort_atoms(const int x1, const int x2, const string& direction = "up");
    
    /** Perform spatial sorting by ordering atoms along Hilbert curve
     *  http://doc.cgal.org/latest/Spatial_sorting/index.html */
    const void sort_spatial();
    
    /** Define the addition of two Mediums */
    Medium& operator +=(Medium &m);

    /** Add data from other Medium to current one */
    const void add(Medium *m);

    /** Reserve memory for data vectors */
    virtual const void reserve(const int n_atoms);

    /** Add atom to the system */
    const void add_atom(const Atom& atom);

    /** Add atom with default id and coordination to the system */
    const void add_atom(const Point3& point);

    /** Write first n_max atoms to file; n_max < 0 writes all the atoms */
    void write(const string &file_name, const int n_max = -1) const;
    
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
    /** Set the marker of i-th atom */
    const void set_marker(const int i, const int m);

    /** Return i-th atom */
    const Atom get_atom(const int i) const;
    /** Return x-, y- and z-coordinate and associated operators for i-th atom */
    const Point3 get_point(const int i) const;
    /** Return x- and y-coordinate and associated operators for i-th atom */
    const Point2 get_point2(const int i) const;
    /** Return ID of i-th atom */
    const int get_id(const int i) const;
    /** Return marker of i-th atom */
    const int get_marker(const int i) const;
    /** Return number of atoms in a Medium */
    const int get_n_atoms() const;

    /** Statistics about system size */
    struct Sizes {
        double xmin;     ///< minimum x-coordinate of atoms
        double xmax;     ///< maximum x-coordinate of atoms
        double ymin;     ///< minimum y-coordinate of atoms
        double ymax;     ///< maximum y-coordinate of atoms
        double zmin;     ///< minimum z-coordinate of atoms
        double zmax;     ///< maximum z-coordinate of atoms
        double zminbox;  ///< minimum z-coordinate of simulation box
        double zmaxbox;  ///< maximum z-coordinate of simulation box
        double xbox;     ///< simulation box size in x-direction
        double ybox;     ///< simulation box size in y-direction
        double zbox;     ///< simulation box size in z-direction
        double xmean;    ///< average value of x-coordinate
        double ymean;    ///< average value of y-coordinate
        double zmean;    ///< average value of z-coordinate
        double xmid;     ///< middle value of x-coordinate
        double ymid;     ///< middle value of y-coordinate
        double zmid;     ///< middle value of z-coordinate
    } sizes;
    vector<Atom> atoms;  ///< vector holding atom coordinates and meta data
protected:


    /** Initialise statistics about the coordinates in Medium */
    const void init_statistics();

    /** Output atom data in .xyz format */
    const void write_xyz(ofstream &outfile, const int n_atoms) const;

    /** Output atom data in .vtk format */
    const void write_vtk(ofstream &outfile, const int n_atoms) const;
    
    /** Get point representation in vtk format */
    virtual const void get_cell_types(ofstream& outfile, const int n_cells) const;

    /** Get data scalar and vector data associated with vtk cells */
    virtual const void get_cell_data(ofstream& outfile, const int n_cells) const;

    /** Get data scalar and vector data associated with vtk nodes */
    virtual const void get_point_data(ofstream& outfile, const int n_atoms) const;

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    virtual const string get_data_string(const int i) const;
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
