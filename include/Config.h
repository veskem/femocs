/*
 * Config.h
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include "Macros.h"

using namespace std;
namespace femocs {

/** Class to read configuration parameters from configuration file */
class Config {
public:

    /** Config constructor */
    Config();

    /** Read the configuration parameters from input script */
    const void read_all(const string& file_name);

    /** Look up the configuration parameter with string argument */
    const int read_command(const string& param, string& arg);

    /** Look up the configuration parameter with boolean argument */
    const int read_command(const string& param, bool& arg);

    /** Look up the configuration parameter with integer argument */
    const int read_command(const string& param, int& arg);

    /** Look up the configuration parameter with double argument */
    const int read_command(const string& param, double& arg);

    /** Look up the configuration parameter with several double arguments */
    const int read_command(const string& param, vector<double>& args);

    /** Print the stored commands and parameters */
    const void print_data();

    string infile;              ///< Path to the file with atom coordinates and types
    string mesh_quality;        ///< Minimum quality (maximum radius-edge ratio) of tetrahedra
    string message;             ///< data string from the host code
    double latconst;            ///< Lattice constant
    double coord_cutoff;        ///< Cut-off distance in Angstroms for Coordination analysis
    int nnn;                    ///< Number of nearest neighbours for given crystal structure
    int nt;                     ///< Number of OpenMP threads
    double neumann;             ///< Value of Neumann boundary condition
    bool postprocess_marking;   ///< Make extra effort to mark correctly the vacuum nodes in shadowed area
    bool refine_apex;           ///< Add elements to the nanotip apex
    double zbox_above;          ///< Space added above the maximum z-coordinate of surface
    double zbox_below;          ///< Space added below the minimum z-coordinate of surface

    /// Distance from surface edge where atoms are picked for rectangularization
    double rmin_rectancularize;

    /** Minimum distance between atoms from current and previous run so that their
     * movement is considered to be sufficiently big to recalculate electric field;
     * 0 turns the check off */
    double distance_tol;

    /// Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
    double radius;

//    /// Factor that is proportional to the extent of surface coarsening; 0 turns coarsening off
//    double coarse_factor;

    /// Factor that is proportional to the extent of surface smoothing; 0 turns smoothing off
    double smooth_factor;

    /// Number of bins in smoother histogram; 1 or less turns off the histogram smoother
    int n_bins;

    /// Factors that are proportional to the extent of surface coarsening; 0 turns corresponding coarsening component off
    struct CoarseFactor {
        double amplitude;
        double r0_cylinder;
        double r0_sphere;
    } cfactor;

private:
    vector<vector<string>> data;  ///< Commands and their parameters found from the input script

    const string comment_symbols = "!#%";
    const string data_symbols = "/*_.0123456789abcdefghijklmnopqrstuvwxyz";

    const void parse_file(const string& file_name);

    /** Initialize configuration parameters */
    const void init_values();

    /** Remove the noise from the beginning of the string */
    const void trim(string& str);
};

} // namespace femocs

#endif /* CONFIG_H_ */
