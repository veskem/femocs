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
    virtual ~Medium() {};

    /** Sort the atoms by their x, y or z coordinate (coord=0|1|2) or radial coordinate(coord=3) */
    void sort_atoms(const int coord, const string& direction="up");

    /** Sort the atoms twice by their x, y or z coordinate */
    void sort_atoms(const int x1, const int x2, const string& direction = "up");
    
    /** Perform spatial sorting by ordering atoms along Hilbert curve
     *  http://doc.cgal.org/latest/Spatial_sorting/index.html */
    void sort_spatial();
    
    /** Append data from other Medium to current one */
    Medium& operator +=(const Medium &m);

    /** Initialize memory for data vectors */
    virtual void reserve(const int n_atoms);
    
    /** Reserve memory for data vectors without erasing already excisting data */
    virtual void resize(const int n_atoms);

    /** Add atom to the system */
    void append(const Atom& atom);

    /** Add atom with default id and coordination to the system */
    void append(const Point3& point);

    /** Write atoms to file. Function is active only when file write is enabled */
    void write(const string &file_name) const;
    
    /** Calculate statistics about the coordinates in Medium */
    void calc_statistics();

    /** Copy statistics from another Medium */
    void copy_statistics(const Medium& m);

    /** Make the atoms coordinates to correspond to another Medium */
    void update_positions(const Medium& medium);

    /** Set the id of i-th atom */
    void set_id(const int i, const int id);
    /** Set the coordinates of i-th atom */
    void set_point(const int i, const Point3& p);
    /** Set the x-coordinate of i-th atom */
    void set_x(const int i, const double x);
    /** Set the y-coordinate of i-th atom */
    void set_y(const int i, const double y);
    /** Set the z-coordinate of i-th atom */
    void set_z(const int i, const double z);
    /** Set the marker of i-th atom */
    void set_marker(const int i, const int m);

    /** Return i-th atom */
    Atom get_atom(const int i) const;
    /** Return x-, y- and z-coordinate and associated operators for i-th atom */
    Point3 get_point(const int i) const;
    /** Return x- and y-coordinate and associated operators for i-th atom */
    Point2 get_point2(const int i) const;
    /** Return ID of i-th atom */
    int get_id(const int i) const;
    /** Return marker of i-th atom */
    int get_marker(const int i) const;
    /** Return number of atoms in a Medium */
    int size() const;

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

        int size() const { return (&zmid)-(&xmin) + 1; }  ///< number of values in Sizes
    } sizes;

    vector<Atom> atoms;  ///< vector holding atom coordinates and meta data
protected:
    vector<array<int,3>> nborbox_indices; ///< neighbour box indices where the point belongs to
    array<int,3> nborbox_size;            ///< # neighbour boxes in x,y,z direction
    vector<int> list;  ///< linked list entries
    vector<int> head;  ///< linked list header

    /**
     * Calculate Verlet neighbour list for atoms by organizing atoms first to linked list.
     * For theory see
     * http://www.acclab.helsinki.fi/~knordlun/moldyn/lecture03.pdf
     * http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
     */
    void calc_verlet_nborlist(vector<vector<int>>& nborlist, const double r_cut, const bool periodic);

    /** Calculate linked list between atoms that holds the information about
     * the region  of simulation cell where the atoms are located.
     * Linked list can be used to calculate efficiently the neighbour list. */
    void calc_linked_list(const double r_cut);

    /** Initialise statistics about the coordinates in Medium */
    void init_statistics();

    /** Get scalar and vector data associated with atoms */
    virtual void get_cell_data(ofstream& outfile) const;

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    virtual string get_data_string(const int i) const;

    /** Get entry to the dat-file; first_line == true gives the header of data */
    virtual string get_global_data(const bool first_line) const;

private:
    /** Output atom data in .xyz format */
    void write_xyz(ofstream &outfile, const int n_atoms) const;

    /** Output atom data in .vtk format */
    void write_vtk(ofstream &outfile, const int n_atoms) const;
    
    /** Output atom data in .ckx format that shows atom coordinates and their types (fixed, surface, bulk etc.)
     * Atom types are the same as in Types struct in Macros.h */
    void write_ckx(ofstream &outfile, const int n_atoms) const;

    /** Append single line of data into dat-file.
     * If the file is empty, the header is written first and then data follows. */
    void write_dat(ofstream &outfile) const;

    void loop_nbor_boxes(vector<vector<int>>& nborlist, const double r_cut2, const int atom);

    void loop_periodic_nbor_boxes(vector<vector<int>>& nborlist, const double r_cut2, const int atom);

    inline int periodic_image(int image, int coordinate) const;
};

} /* namespace femocs */

#endif /* MEDIUM_H_ */
