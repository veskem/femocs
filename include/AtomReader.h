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
#include "Config.h"

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
    void generate_nanotip(double height, double radius, double latconst);

    /** Import atom coordinates and types from a file and check their rmsd
     * @param file_name  path to input file with atomic data
     * @param add_noise  add random noise to the imported atom coordinates to emulate real simulation
     */
    bool import_file(const string &file_name, const bool add_noise=false);

    /** Import atoms coordinates from PARCAS and check their rmsd
     * @param n_atoms     number of imported atoms
     * @param coordinates vector of atomistic coordinates in PARCAS units; x0=c[0], y0=c[1], z[0]=c[2], x1=c[3], etc
     * @param box         vector of MD simulation box sizes in Angstroms; box[0]->x_box, box[1]->y_box, box[2]->z_box
     */
    bool import_parcas(const int n_atoms, const double* coordinates, const double* box);

    /** Import atom coordinates and types and check their rmsd */
    bool import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types);

    /** Calculate coordination for all the atoms by using PARCAS neighbour list
     * or, if it's missing, building first the Verlet neighbour list */
    void calc_coordinations(const int* parcas_nborlist=NULL);

    /** Calculate coordination for all the atoms by using PARCAS neighbour list
     * or, if it's missing, building first the Verlet neighbour list.
     * Before doing so, update nnn, lattice constant and coordination cut-off radius
     * by calculating radial distribution function. */
    void calc_rdf_coordinations(const int* parcas_nborlist=NULL);

    /** Calculate pseudo-coordination for all the atoms using the atom types */
    void calc_pseudo_coordinations();

    /** Rebuild list of close neighbours and run cluster analysis.
     * Atoms are grouped into clusters using density-based spatial clustering technique
     * http://codereview.stackexchange.com/questions/23966/density-based-clustering-of-image-keypoints
     * https://en.wikipedia.org/wiki/DBSCAN */
    void calc_clusters(const int* parcas_nborlist=NULL);

    /** Extract atom types from calculated atom coordinations */
    void extract_types();

    /** Store the atom coordinates from current run */
    void save_current_run_points(const double eps);

    /** Store commonly used data */
    void store_data(const Config& conf);

    /** Obtain rms-distance the atoms have moved since the last full iteration.
     * BIG values (>10^308) indicate that current and previous iterations are not comparable. */
    double get_rmsd() const { return data.rms_distance; }

    /** Return number of atom that are detached from the big system */
    int get_n_detached() const { return data.n_detached; }

    /** Return factors to covert SI units to Parcas ones */
    Vec3 get_si2parcas_box() const;

    /** Return factors to covert Parcas units to SI ones */
    Vec3 get_parcas2si_box() const;

    /** Obtain the printable statistics */
    friend ostream& operator <<(ostream &os, const AtomReader& ar) {
        os << fixed << setprecision(3)
                << "nnn=" << ar.data.nnn
                << ", latconst=" << ar.data.latconst
                << ", coord_cutoff=" << ar.data.coord_cutoff
                << ", cluster_cutoff=" << ar.data.cluster_cutoff
                << ", #evap|clust atoms=" << ar.data.n_evaporated << "|" << ar.data.n_detached;
        return os;
    }

private:
    vector<int> cluster;            ///< id of cluster the atom is located
    vector<int> coordination;       ///< coordinations of atoms
    vector<int> previous_types;     ///< atom types from previous run
    vector<Point3> previous_points; ///< atom coordinates from previous run
    vector<vector<int>> nborlist;   ///< list of closest neighbours
    Vec3 simubox;                   ///< MD simulation box dimensions; needed to convert SI units to Parcas one

    struct Data {
        double distance_tol=0;        ///< min rms distance atoms must move to consider their movement significantly big
        double rms_distance=0;        ///< rms distance between atoms from previous and current run
        double cluster_cutoff=0;      ///< cluster analysis cut-off radius
        double coord_cutoff=0;        ///< coordination analysis cut-off radius
        double latconst=0;            ///< lattice constant
        unsigned int nnn=0;           ///< effective number of nearest neighbours
        unsigned int n_detached=0;    ///< number of atoms that are detached from the big structure
        unsigned int n_evaporated=0;  ///< number of atoms that are evaporated from the big structure
    } data;

    /** Import atoms from different types of file.
     * @param file_name - path to file with atomic data
     */
    void import_xyz(const string& file_name);
    void import_ckx(const string& file_name);

    /** Reserve memory for data vectors */
    void reserve(const int n_atoms);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;

    /** Calculate list of close neighbours using Parcas diagonal neighbour list */
    void calc_nborlist(const double r_cut, const int* parcas_nborlist);

    /** Calculate list of close neighbours using already existing list with >= cut-off radius */
    void recalc_nborlist(const double r_cut);

    /** Calculate the radial distribution function (rdf) in a periodic isotropic system.
     *  Source of inspiration: https://github.com/anyuzx/rdf
     *  Author: Guang Shi, Mihkel Veske
    */
    void calc_rdf(const int n_bins, const double r_cut);
    void calc_rdf_peaks(vector<double>& peaks, const vector<double>& rdf, const double bin_width);

    /** Calculate the root mean square average distance the atoms have moved
     * between previous and current run */
    bool calc_rms_distance();
};

} /* namespace femocs */

#endif /* ATOMREADER_H_ */
