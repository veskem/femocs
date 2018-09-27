/*
 * Config.h
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include "Macros.h"
#include <set>

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

    /** Look up the configuration parameter with unsigned integer argument */
    int read_command(string param, unsigned int& arg);

    /** Look up the configuration parameter with integer argument */
    int read_command(string param, int& arg);

    /** Look up the configuration parameter with double argument */
    int read_command(string param, double& arg);

    /** Look up the configuration parameter with several double arguments */
    int read_command(string param, vector<double>& args);

    /** Print the stored commands and parameters */
    void print_data();

    /** Various paths */
    struct Path {
        string extended_atoms;      ///< Path to the file with atoms forming the extended surface
        string infile;              ///< Path to the file with atom coordinates and types
        string mesh_file;           ///< Path to the triangular and tetrahedral mesh data
    } path;

    /** User specific preferences */
    struct Behaviour {
        string verbosity;           ///< Verbose mode: mute, silent, verbose
        string project;             ///< Type of project to be called
        int n_writefile;            ///< Number of time steps between writing output files; 0 turns writing off
        int interpolation_rank;     ///< Rank of the solution interpolation; 1-linear tetrahedral, 2-quadratic tetrahedral, 3-linear hexahedral
        double write_period;        ///< Write files every write_period (in fs)
        double timestep_fs;         ///< Total time evolution within a FEMOCS run call [fs]
        double mass;                ///< Atom mass [amu]
        unsigned int rnd_seed;      ///< Seed for random number generator
        unsigned int n_omp_threads; ///< Number of opened OpenMP threads
    } behaviour;

    /** Enable or disable various support features */
    struct Run {
        bool cluster_anal;          ///< Enable cluster analysis
        bool apex_refiner;          ///< Add elements to the nanotip apex
        bool rdf;                   ///< Re-calculate lattice constant and coordination analysis parameters using radial distribution function
        bool output_cleaner;        ///< Clear output folder before the run
        bool surface_cleaner;       ///< Clean surface by measuring the atom distance from the triangular surface
        bool field_smoother;        ///< Replace nodal field with the average of its neighbouring nodal fields
    } run;

    /** Sizes related to mesh, atoms and simubox */
    struct Geometry {
        string mesh_quality;        ///< Minimum quality (maximum radius-edge ratio) of tetrahedra
        string element_volume;      ///< Maximum volume of tetrahedra
        int nnn;                    ///< Number of nearest neighbours for given crystal structure
        double latconst;            ///< Lattice constant
        double coordination_cutoff; ///< Cut-off distance for coordination analysis [same unit as latconst]
        double cluster_cutoff;      ///< Cut-off distance for cluster analysis [same unit as latconst]; if 0, cluster analysis uses coordination_cutoff instead
        double charge_cutoff;       ///< Cut-off distance for calculating Coulomb forces [same unit as latconst]
        double surface_thickness;   ///< Maximum distance the surface atom is allowed to be from surface mesh [same unit as latconst]; 0 turns check off
        double box_width;           ///< Minimal simulation box width [tip height]
        double box_height;          ///< Simulation box height [tip height]
        double bulk_height;         ///< Bulk substrate height [lattice constant]

        /// Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
        double radius;
        double height;              ///< height of generated artificial nanotip in the units of radius
    } geometry;

    /** All kind of tolerances */
    struct Tolerance {
        double charge_min; ///< Min ratio face charges are allowed to deviate from the total charge
        double charge_max; ///< Max ratio face charges are allowed to deviate from the total charge
        double field_min;  ///< Min ratio numerical field can deviate from analytical one
        double field_max;  ///< Max ratio numerical field can deviate from analytical one

        /** Minimum rms distance between atoms from current and previous run so that their
         * movement is considered to be sufficiently big to recalculate electric field;
         * 0 turns the check off */
        double distance;
    } tolerance;

    /** Parameters for solving field equation */
    struct Field {
        double E0;              ///< Value of long range electric field (Active in case of Neumann anodeBC
        double ssor_param;      ///< Parameter for SSOR preconditioner in DealII
        double phi_error;       ///< Maximum allowed electric potential error
        int n_phi;              ///< Maximum number of Conjugate Gradient iterations in phi calculation
        double V0;              ///< Applied voltage at the anode (active in case of SC emission and Dirichlet anodeBC
        string anode_BC;         ///< Type of anode boundary condition (Dirichlet or Neumann)
        string solver;          ///< Type of field equation to be solved; laplace or poisson
        int element_degree;     ///< Degree of Finite elements (1: linear, 2: quadratic, 3: cubic ...
    } field;

    /** Heating module configuration parameters */
    struct Heating {
        string mode;                ///< Method to calculate current density and temperature; none, stationary or transient
        string rhofile;             ///< Path to the file with resistivity table
        double lorentz;             ///< Lorentz number (Wiedemenn-Franz law)
        double t_ambient;           ///< Ambient temperature in heat calculations
        double t_error;             ///< Maximum allowed temperature error in Newton iterations
        int n_newton;               ///< Maximum number of Newton iterations
        int n_cg;                   ///< Max # Conjugate-Gradient iterations
        double cg_tolerance;        ///< Solution accuracy in Conjugate-Gradient solver
        double ssor_param;          ///< Parameter for SSOR preconditioner in DealII. Its fine tuning optimises calculation time.
        double delta_time;          ///< Timestep of time domain integration [sec]
        double dt_max;              ///< Maximum allowed timestep for heat convergence run
        double tau;                 ///< Time constant in Berendsen thermostat
        string assemble_method;     ///< Method to assemble system matrix for solving heat equation; euler or crank_nicolson
    } heating;

    struct Emission {
        double work_function;       ///< Work function [eV]
        bool blunt;                 ///< Force blunt emitter approximation (good for big systems)
        bool cold;                  ///< force cold field emission approximation (good for low temperatures)
        double omega_SC;            ///< Voltage correction factor for SC calculation (negative for ignoring SC)
        double SC_error;            ///< convergence criterion for SC error
        double Vappl_SC;            ///< Applied voltage used for SC calculations (overrides Vappl * omega_SC)
        string SC_mode;             ///< Mode by which the SC is taken into account. global or local
    } emission;

    /** Parameters related to atomic force calculations */
    struct Force {
        string mode;                ///< Forces to be calculated; lorentz, all, none
    } force;

    /** Smooth factors for surface faces, surface atoms and charges */
    struct Smoothing {
        string algorithm;    ///< surface mesh smoother algorithm; none, laplace or fujiwara
        int n_steps;         ///< number of surface mesh smoothing iterations
        double lambda_mesh;  ///< lambda parameter in surface mesh smoother
        double mu_mesh;      ///< mu parameter in surface mesh smoother
        double beta_atoms;   ///< extent of surface smoothing; 0 turns smoothing off
        double beta_charge;  ///< extent of charge smoothing; 0 turns smoothing off
    } smoothing;

    /** Factors that are proportional to the extent of surface coarsening; 0 turns corresponding coarsening component off */
    struct CoarseFactor {
        double amplitude;
        int r0_cylinder;
        int r0_sphere;
        double exponential;
    } cfactor;

    /** Particle In Cell module configuration */
    struct PIC {
        string mode;      ///< Pic mode (transient, converge or none)

        /** Maximum PIC timestep [fs].
         * The actual PIC timestep will be smaller
         * such that it is an integer fraction of the MD timestep. */
        double dt_max;
        double Wsp_el;        ///< Superparticle weight for electrons
        bool fractional_push; ///< Do fractional timestep push when injecting electrons?
        bool coll_coulomb_ee; ///< Do 2e->2e Coulomb collisions?
    } pic;
    
    /** Parameters related to SpaceCharge project */
    struct SpaceCharge {
        vector<double> apply_factors; ///< run for multiple applied E0 (or V0) multiplied by the factors
        double convergence;     ///< Relative error in current for convergence criterion.
                                ///< It is compared with the corresponding std in the pic convergence step.
        vector<double> I_pic;   ///< Current target for finding Applied SC voltage (for SC calculations)
    } SC;

private:
    vector<vector<string>> data;          ///< commands and their arguments found from the input script

    const string comment_symbols = "!#%";
    const string data_symbols = "+-/*_.0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ()";

    /** Check for the obsolete commands from the buffered commands */
    void check_obsolete(const string& file_name);

    /** Check for the obsolete commands that are similar to valid ones */
    void check_obsolete(const string& command, const string& substitute);

    /** Read the commands and their arguments from the file and store them into the buffer */
    void parse_file(const string& file_name);

    /** Remove the noise from the beginning of the string */
    void trim(string& str);
};

} // namespace femocs

#endif /* CONFIG_H_ */
