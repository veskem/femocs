/*
 * AtomReader.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include "Macros.h"
#include "Medium.h"

using namespace std;
namespace femocs {

/** Class to import atoms from atomistic simulation and to divide them into different categories */
class AtomReader: public Medium {
public:
    /** Constructor for AtomReader */
    AtomReader();

    void generate_nanotip(const double height, const double radius, const double latconst);

    /**
     * Function to import file with atom coordinates and types
     * @param file_name  path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     */
    void import_file(const string &file_name);

    /** Function to transform atomic data from Parcas format into AtomReader one */
    void import_parcas(const int n_atoms, const double* coordinates, const double* box);

    /** Function to transform atomic data from Helmod format into AtomReader one */
    void import_helmod(const int n_atoms, const double* x, const double* y, const double* z, const int* types);

    /** Function to transform atomic data from Kimocs format into AtomReader one */
    void import_kimocs();

    /** Calculate coordination for all the atoms using Parcas neighbour list */
    void calc_coordinations(const int nnn, const double r_cut, const int* parcas_nborlist);

    /** Calculate coordination for all the atoms using brute force */
    void calc_coordinations(const int nnn, const double r_cut);

    /** Calculate pseudo-coordination for all the atoms using the atom types */
    void calc_coordinations(const int nnn);

    /** Group atoms into clusters using density-based spatial clustering technique
     * http://codereview.stackexchange.com/questions/23966/density-based-clustering-of-image-keypoints
     * https://en.wikipedia.org/wiki/DBSCAN */
    void calc_clusters();

    /** Rebuild list of close neighbours using brute force technique and run cluster analysis */
    void calc_clusters(const int nnn, const double r_cut);

    /** Rebuild list of close neighbours using Parcas neighbourlist and run cluster analysis */
    void calc_clusters(const int nnn, const double r_cut, const int* parcas_nborlist);

    /** Function to detect evaporated atoms by their coordinations;
     * an atom is considered evaporated if it's coordination is between the defined limits. */
    void check_coordinations();

    /** Calculate the number of atoms in clusters
     * @param print  if clusters detected, print the number atoms in clusters to the console
     */
    int check_clusters(const bool print);

    /** Extract atom types from calculated atom coordinations */
    void extract_types(const int nnn, const double latconst);

    /** Redefine the min and max values for z-coordinates */
    void resize_box(const double zmin, const double zmax);
    
    /** Redefine the min and max values for x, y and z - coordinates */
    void resize_box(const double xmin, const double xmax, const double ymin, const double ymax,
            const double zmin, const double zmax);

    /** Calculate the root mean square average distance the atoms have moved
     * between previous and current run */
    double calc_rms_distance(const double eps);
    
    /** Store the atom coordinates from current run */
    void save_current_run_points(const double eps);

    /** Get the closest neighbours of i-th atom */
    const vector<int>& get_neighbours(const int i) const;

    /** Return the size of neighbour list */
    int get_nborlist_size() const;

    double rms_distance;            ///< rms distance between atoms from previous and current run

private:
    vector<int> cluster;            ///< id of cluster the atom is located
    vector<int> coordination;       ///< coordinations of atoms
    vector<int> previous_types;     ///< atom types from previous run
    vector<Point3> previous_points; ///< atom coordinates from previous run
    vector<vector<int>> nborlist;   ///< list of closest neighbours

    /**
     * Functions to import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    void import_xyz(const string& file_name);
    void import_ckx(const string& file_name);
    void import_dump(const string& file_name);

    /** Reserve memory for data vectors */
    void reserve(const int n_atoms);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;

    /** Calculate list of close neighbours using Parcas diagonal neighbour list */
    void calc_nborlist(const int nnn, const double r_cut, const int* parcas_nborlist);

    /** Calculate list of close neighbours using brute force technique */
    void calc_nborlist(const int nnn, const double r_cut);
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
