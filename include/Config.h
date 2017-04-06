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

/** Class to initialize and read configuration parameters from configuration file */
class Config {
public:

    /** Config constructor initializes configuration parameters */
    Config();

    /** Read the configuration parameters from input script */
    void read_all(const string& file_name);

    /** Look up the configuration parameter with string argument */
    int read_command(string param, string& arg);

    /** Look up the configuration parameter with boolean argument */
    int read_command(string param, bool& arg);

    /** Look up the configuration parameter with integer argument */
    int read_command(string param, int& arg);

    /** Look up the configuration parameter with double argument */
    int read_command(string param, double& arg);

    /** Look up the configuration parameter with several double arguments */
    int read_command(string param, vector<double>& args);

    /** Print the stored commands and parameters */
    void print_data();

    string extended_atoms;      ///< Path to the file with atoms forming the extended surface
    string atom_file;           ///< Path to the file with atom coordinates and types
    string mesh_quality;        ///< Minimum quality (maximum radius-edge ratio) of tetrahedra
    string element_volume;      ///< Maximum volume of tetrahedra
    string message;             ///< Data string from the host code
    double latconst;            ///< Lattice constant
    double coord_cutoff;        ///< Cut-off distance in Angstroms for Coordination analysis
    double surface_thichness;   ///< Maximum distance the surface atom is allowed to be from surface mesh [angstrom]; 0 turns check off
    int nnn;                    ///< Number of nearest neighbours for given crystal structure
    int nt;                     ///< Number of OpenMP threads
    double neumann;             ///< Value of Neumann boundary condition
    double E0;                  ///< Value of long range electric field
    bool cluster_anal;          ///< Enable cluster analysis
    bool postprocess_marking;   ///< Make extra effort to mark correctly the vacuum nodes in shadowed area
    bool refine_apex;           ///< Add elements to the nanotip apex
    bool heating;               ///< Turn ON 3D current density and temperature calculations
    double box_width;           ///< Minimal simulation box width [tip height]
    double box_height;          ///< Simulation box height [tip height]
    double bulk_height;         ///< Bulk substrate height [lattice constant]
    double t_error;             ///< Maximum allowed temperature error in Newton iterations
    int n_newton;               ///< Maximum number of Newton iterations
    double ssor_param;          ///< Parameter for SSOR preconditioner in DealII
    double phi_error;           ///< Maximum allowed electric potential error
    int n_phi;                  ///< Maximum number of Conjugate Gradient iterations in phi calculation
    bool clear_output;          ///< Clear output folder before the run
    bool use_histclean;         ///< Clean the solution with histogram cleaner
    int n_writefile;            ///< Number of time steps between writing output files; 0 turns writing off
    string verbose_mode;        ///< Verbose mode: mute, silent, verbose
    double force_factor;        ///< Factor determining the relationship between force and charge*elfield

    /** Method to clean the surface atoms
     * voronois - use Voronoi cells
     * faces - measure distance from surface faces
     * none - do not use the cleaner (because it is guaranteed to be clean)
     */
    string surface_cleaner;

    /** Minimum distance between atoms from current and previous run so that their
     * movement is considered to be sufficiently big to recalculate electric field;
     * 0 turns the check off */
    double distance_tol;

    /// Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
    double radius;

    /// Factor that is proportional to the extent of surface smoothing; 0 turns smoothing off
    double surface_smooth_factor;

    /// Factor that is proportional to the extent of charge smoothing; 0 turns smoothing off
    double charge_smooth_factor;

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

    void parse_file(const string& file_name);

    /** Remove the noise from the beginning of the string */
    void trim(string& str);
};

} // namespace femocs

#endif /* CONFIG_H_ */
