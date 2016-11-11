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

/* This is a simple reader; more complicated ones could be obtained by
 * using dedicated library, for example libconfig,
 * http://www.hyperrealm.com/libconfig/
 *
 * See more
 * http://stackoverflow.com/questions/6892754/creating-a-simple-configuration-file-and-parser-in-c
 * */
class Config {
public:
    Config();

    const void read_all(const string& file_name);
    const void read_parameter(const string& param, string& value);
    const void read_parameter(const string& param, bool& value);
    const void read_parameter(const string& param, int& value);
    const void read_parameter(const string& param, double& value);

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
    double significant_distance;

    /// Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
    double radius;

    /// Factor that is proportional to the extent of surface coarsening; 0 turns coarsening off
    double coarse_factor;

    /// Factor that is proportional to the extent of surface smoothing; 0 turns smoothing off
    double smooth_factor;

    /// Number of bins in smoother histogram; 1 or less turns off the histogram smoother
    int n_bins;
private:
    string infile;              ///< Path to input script
    const void init_values();
    const void trim(string& str);
};

} // namespace femocs

#endif /* CONFIG_H_ */
