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
#include "Surface.h"

using namespace std;
namespace femocs {

/** Class to import atoms from atomistic simulation and to divide them into different categories */
class AtomReader: public Medium {
public:
    /** Constructor for AtomReader */
    AtomReader();

    /** Extract atom with desired types and clean it from lonely atoms;
     * atom is considered lonely if its coordination is lower than threshold.
     * @param surface  Surface where the extracted atoms will be written
     * @param type     type of the atoms that will be read
     * @param invert   if true all the atoms except the 'type'-ones will be stored
     */
    void extract(Surface& surface, const int type, const bool invert=false);

    /** Generate nanotip with high rotational symmetry and without crystallographic faceting */
    void generate_nanotip(const double height, const double radius, const double latconst);

    /**
     * Function to import file with atom coordinates and types
     * @param file_name  path to input file with atomic data in .xyz (PARCAS), .dump (LAMMPS) or .ckx (KIMOCS) format
     * @param add_noise  add random noise to the imported atom coordinates to emulate real simulation
     */
    void import_file(const string &file_name, const bool add_noise=false);

    /** Function to transform atomic data from Parcas format into AtomReader one */
    void import_parcas(const int n_atoms, const double* coordinates, const double* box);

    /** Function to transform atomic data from Helmod format into AtomReader one */
    void import_helmod(const int n_atoms, const double* x, const double* y, const double* z, const int* types);

    /** Function to transform atomic data from Kimocs format into AtomReader one */
    void import_kimocs();

    /** Calculate coordination for all the atoms using Parcas neighbour list */
    void calc_coordinations(const int nnn, const double coord_cutoff, const int* parcas_nborlist);

    /** Calculate coordination for all the atoms using neighbour list.
     * Before doing so, update nnn, lattice constant and coordination cut-off radius
     * by calculating radial distribution function */
    void calc_coordinations(int& nnn, double& latconst, double& coord_cutoff, const int* parcas_nborlist);

    /** Calculate coordination for all the atoms using brute force */
    void calc_coordinations(const int nnn, const double coord_cutoff);
    void calc_coordinations(int& nnn, double& latconst, double& coord_cutoff);

    /** Calculate pseudo-coordination for all the atoms using the atom types */
    void calc_coordinations(const int nnn);

    /** Rebuild list of close neighbours using brute force technique and run cluster analysis */
    void calc_clusters(const double cluster_cutoff, const double coord_cutoff);

    /** Rebuild list of close neighbours using Parcas neighbourlist and run cluster analysis */
    void calc_clusters(const int nnn, const double cluster_cutoff, const double coord_cutoff, const int* parcas_nborlist);

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

    /** Return factors to covert SI units to Parcas ones */
    Vec3 get_si2parcas_box() const;

    /** Return factors to covert Parcas units to SI ones */
    Vec3 get_parcas2si_box() const;

    double rms_distance;            ///< rms distance between atoms from previous and current run

private:
    vector<int> cluster;            ///< id of cluster the atom is located
    vector<int> coordination;       ///< coordinations of atoms
    vector<int> previous_types;     ///< atom types from previous run
    vector<Point3> previous_points; ///< atom coordinates from previous run
    vector<vector<int>> nborlist;   ///< list of closest neighbours
    Vec3 simubox;                   ///< MD simulation box dimensions; needed to convert SI units to Parcas one

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

    /** Calculate list of close neighbours using already existing list with >= cut-off radius */
    void recalc_nborlist(const double r_cut);

    /** Function to calculate the radial distribution function in a periodic condition for isotropic system.
     *  Source of inspiration: https://github.com/anyuzx/rdf
     *  Author: Guang Shi, Mihkel Veske
    */
    void calc_rdf(int & nnn, double& latconst, double& coord_cutoff, const int n_bins, const double r_cut);
    void calc_rdf_peaks(vector<double>& peaks, const vector<double>& rdf, const double bin_width);

    /** Group atoms into clusters using density-based spatial clustering technique
     * http://codereview.stackexchange.com/questions/23966/density-based-clustering-of-image-keypoints
     * https://en.wikipedia.org/wiki/DBSCAN */
    void calc_clusters();

    /** Using the previously calculated neighbour list, calculate the coordinations of atoms */
    void calc_coordinations();
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
