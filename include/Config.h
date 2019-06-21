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

    /** Constructor initializes configuration parameters */
    Config();

    /** Read the configuration parameters from input script */
    void read_all(const string& file_name, bool full_run=true);

    /** Using the stored path, read the configuration parameters from input script */
    void read_all();

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
        string restart_file;        ///< Path to the restart file
    } path;

    /** User specific preferences */
    struct Behaviour {
        string verbosity;           ///< Verbose mode: mute, silent, verbose
        string project;             ///< Type of project to be called; runaway, ...
        int n_restart;              ///< # time steps between writing restart file; 0 disables write
        int restart_multiplier;     ///< After n_write_restart * restart_multiplier time steps, restart file will be copied to separate file
        int n_write_log;            ///< # time steps between writing log file; <0: only last time step, 0: no write, >0: only every n-th
        int n_read_conf;            ///< # time steps between re-reading configuration values from file; 0 turns re-reading off
        int interpolation_rank;     ///< Rank of the solution interpolation; 1-linear tetrahedral, 2-quadratic tetrahedral, 3-linear hexahedral
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
        bool smooth_updater;        ///< Force and field will be evaluated, fully or partially, every time step
    } run;

    /** Sizes related to mesh, atoms and simubox */
    struct Geometry {
        int nnn;                    ///< Number of nearest neighbours for given crystal structure
        double latconst;            ///< Lattice constant
        double coordination_cutoff; ///< Cut-off distance for coordination analysis [same unit as latconst]
        double cluster_cutoff;      ///< Cut-off distance for cluster analysis [same unit as latconst]; if 0, cluster analysis uses coordination_cutoff instead
        double charge_cutoff;       ///< Cut-off distance for calculating Coulomb forces [same unit as latconst]
        double surface_thickness;   ///< Maximum distance the surface atom is allowed to be from surface mesh [same unit as latconst]; 0 turns check off
        double box_width;           ///< Minimal simulation box width [tip height]
        double box_height;          ///< Simulation box height [tip height]
        double bulk_height;         ///< Bulk substrate height [lattice constant]
        double radius;              ///< Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
        double theta;               ///< Apex angle of coarsening cone [degree]
        double height;              ///< height of generated artificial nanotip in the units of radius
        /** Minimum rms distance between atoms from current and previous run so that their
         * movement is considered to be sufficiently big to recalculate electric field;
         * 0 turns the check off */
        double distance_tol;
        double beta_atoms;          ///< Extent of surface smoothing; 0 turns smoothing off
    } geometry;

    /** All kind of tolerances */
    struct Tolerance {
        double charge_min; ///< Min ratio face charges are allowed to deviate from the total charge
        double charge_max; ///< Max ratio face charges are allowed to deviate from the total charge
        double field_min;  ///< Min ratio numerical field can deviate from analytical one
        double field_max;  ///< Max ratio numerical field can deviate from analytical one
    } tolerance;

    /** Parameters for solving field equation */
    struct Field {
        double E0;             ///< Value of long range electric field (Active in case of Neumann anodeBC
        double ssor_param;     ///< Parameter for SSOR preconditioner in DealII
        double cg_tolerance;   ///< Maximum allowed electric potential error
        int n_cg;              ///< Maximum number of Conjugate Gradient iterations in phi calculation
        double V0;             ///< Applied voltage at the anode (active in case of SC emission and Dirichlet anodeBC
        string anode_BC;       ///< Type of anode boundary condition (Dirichlet or Neumann)
        string mode;           ///< Mode to run field solver; laplace, transient or converge; transient or converge activate PIC solver
        double V_min;          ///< Minimum allowed electric potential [V]
        double V_max;          ///< Maximum allowed electric potential [V]
    } field;

    /** Heating module configuration parameters */
    struct Heating {
        string mode;                ///< Method to calculate current density and temperature; none, stationary or transient
        string rhofile;             ///< Path to the file with resistivity table
        double lorentz;             ///< Lorentz number (Wiedemenn-Franz law)
        double t_ambient;           ///< Ambient temperature in heat calculations
        int n_cg;                   ///< Max # Conjugate-Gradient iterations
        double cg_tolerance;        ///< Solution accuracy in Conjugate-Gradient solver
        double ssor_param;          ///< Parameter for SSOR preconditioner in DealII. Its fine tuning optimises calculation time.
        double delta_time;          ///< Timestep of time domain integration [fs]
        double dt_max;              ///< Maximum allowed timestep for heat convergence run
        double tau;                 ///< Time constant in Berendsen thermostat
        double T_min;               ///< Minimum allowed temperature [K]
        double T_max;               ///< Maximum allowed temperature [K]
    } heating;

    /** Field emission module parameters */
    struct Emission {
        double work_function; ///< Work function [eV]
        bool blunt;           ///< Force blunt emitter approximation (good for big systems)
        bool cold;            ///< force cold field emission approximation (good for low temperatures)
        double omega;         ///< Voltage correction factor for SC-limited emission calculation; <= 0 ignores SC
        double J_min;         ///< Minimum current density from single face [amps/Angstrom^2]
        double J_max;         ///< Maximum current density from single face [amps/Angstrom^2]
    } emission;

    /** Parameters related to atomic force calculations */
    struct Force {
        string mode;         ///< forces to be calculated; lorentz, all, none
        double beta;         ///< extent of charge smoothing; 0 turns smoothing off
    } force;

    /** Mesh geometry data */
    struct Mesh {
        string quality;      ///< minimum quality (maximum radius-edge ratio) of tetrahedra
        string volume;       ///< maximum volume of tetrahedra
        string algorithm;    ///< surface mesh smoother algorithm; none, laplace or fujiwara
        int n_steps;         ///< number of surface mesh mesh iterations
        double lambda;       ///< lambda parameter in surface mesh smoother
        double mu;           ///< mu parameter in surface mesh smoother
        double coplanarity;  ///< minimum allowed tetrahedron coplanarity after mesh smoothing
    } mesh;

    /** Factors that are proportional to the extent of surface coarsening; 0 turns corresponding coarsening component off */
    struct CoarseFactor {
        double amplitude;   ///< coarsening factor outside the warm region
        int r0_cylinder;    ///< minimum distance between atoms in nanotip outside the apex
        int r0_sphere;      ///< minimum distance between atoms in nanotip apex
        double exponential; ///< coarsening rate; min distance between coarsened atoms outside the warm region is d_min ~ pow(|r1-r2|, exponential)
    } cfactor;

    /** Particle In Cell module configuration */
    struct PIC {
        /** Maximum PIC timestep [fs]. The actual PIC timestep will be smaller,
         * so that it is an integer fraction of the MD timestep. */
        double dt_max;
        double weight_el;      ///< Electron superparticle weight
        bool fractional_push;  ///< Do fractional timestep push when injecting electrons?
        bool collide_ee;       ///< Do 2e->2e Coulomb collisions?
        bool periodic;         ///< SP-s will be mapped back to simubox in x,y-direction?
        double landau_log;     ///< Landau logarithm
        unsigned int max_injected; ///< Max nr of super particles injected during one step
    } pic;
    
    /** Parameters related to SpaceCharge project */
    struct SpaceCharge {
        /** Relative error in current for convergence criterion.
         * It is compared with the corresponding std in the pic convergence step. */
        double convergence;
        vector<double> apply_factors; ///< run for multiple applied E0 (or V0) multiplied by the factors
        vector<double> I_pic;   ///< Current target for finding Applied SC voltage (for SC calculations)
    } scharge;

private:
    vector<vector<string>> data; ///< commands and their arguments found from the input script
    string file_name;            ///< path to configuration file

    const string comment_symbols = "!#%";
    const string data_symbols = "+-/*_.0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ()";

    /** Check for the obsolete commands from the buffered commands */
    void check_obsolete(const string& file_name);

    /** Check for the obsolete commands that are similar to valid ones */
    void check_changed(const string& command, const string& substitute);

    /** Read the commands and their arguments from the file and store them into the buffer */
    void parse_file(const string& file_name);

    /** Remove the noise from the beginning of the string */
    void trim(string& str);
};

} // namespace femocs

#endif /* CONFIG_H_ */
