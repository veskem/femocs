/*
 * ProjectRunaway.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef PROJECTRUNAWAY_H_
#define PROJECTRUNAWAY_H_

#include "GeneralProject.h"
#include "Surface.h"
#include "SolutionReader.h"
#include "EmissionReader.h"
#include "Interpolator.h"
#include "PhysicalQuantities.h"
#include "PoissonSolver.h"
#include "CurrentHeatSolver.h"
#include "Pic.h"


using namespace std;
namespace femocs {

/**
 * Class for calculating electric field and heating effects
 * around a nanostructure that fits into central cylinder
 * by solving Laplace or Poisson equation
 * and taking into account change of temperature due to field emission.
 */
class ProjectRunaway : public GeneralProject {
public:
    ProjectRunaway(AtomReader &reader, Config &conf);
    ~ProjectRunaway() {}

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const int timestep=-1, const double time=-1);

    /** Export the solution on the imported atomistic points.
     * @param n_points   # of first imported atoms where data exported; 0 disables the export
     * @param data_type  label of the data to be exported. Lower-case label invokes data appending, while upper-case letters clear the previous data.
     * @param data       array where results are written. Vector data is exported coordinate-wise, i.e in a form x1,y1,z1,x2,y2,...
     */
    int export_data(double* data, const int n_points, const string &data_type);

    /** Export the integer data on the imported atomistic points.
     * @param n_points   # of first imported atoms where data exported; 0 disables the export
     * @param data_type  label of the data to be exported. Lower-case label invokes data appending, while upper-case letters clear the previous data.
     * @param data       array where results are written.
     */
    int export_data(int* data, const int n_points, const string& data_type);

    /** Export the solution on the specified points.
      * @param n_points      # of first imported atoms where data exported; 0 disables the export
      * @param data_type     label of the data to be exported. See Labels class for a list of possible cmd-s
      * @param near_surface  data is located on or near the surface.
      * @param x,y,z         coordinates of the points where interpolation is performed
      * @param data          array where results are written. Vector data is exported coordinate-wise, i.e in a form x1,y1,z1,x2,y2,...
      * @param flag          array showing whether specified point was located inside the mesh (1) or not (0)
      */
    int interpolate(double* data, int* flag,
            const int n_points, const string &data_type, const bool near_surface,
            const double* x, const double* y, const double* z);

    /** Read and generate simulation data to continue running interrupted simulation */
    int restart(const string &path_to_file);

protected:
    bool fail;                  ///< If some process failed
    double t0;                  ///< CPU timer
    bool mesh_changed;          ///< True if new mesh has been created
    bool first_run;             ///< True only as long as there is no full run
    double last_heat_time;      ///< Last time heat was updated
    double last_pic_time;       ///< Last time PIC solver was called
    int last_restart_ts;        ///< Last time step reset file was written
    int restart_cntr;           ///< How many restart files have been written

    Interpolator vacuum_interpolator;  ///< data & operations for interpolating field & potential in vacuum
    Interpolator bulk_interpolator;    ///< data & operations for interpolating current density & temperature in bulk

    Surface dense_surf;       ///< non-coarsened surface atoms
    Surface extended_surf;    ///< atoms added for the surface atoms

    FieldReader fields;       ///< fields & potentials on surface atoms
    HeatReader  temperatures; ///< temperatures & current densities on bulk atoms
    ForceReader forces;       ///< forces & charges on surface atoms

    FieldReader surface_fields;       ///< fields on surface hex face centroids
    HeatReader  surface_temperatures; ///< temperatures & current densities on surface hex face centroids
    HeatReader  heat_transfer;        ///< temperatures on new mesh dofs interpolated on old solution space

    PhysicalQuantities phys_quantities; ///< quantities used in heat calculations
    PoissonSolver<3> poisson_solver;    ///< Poisson equation solver
    CurrentHeatSolver<3> ch_solver;     ///< transient currents and heating solver

    EmissionReader emission;          ///< emission data on centroids of surface quadrangles
    Pic<3> pic_solver;                       ///< class for solving Poisson equation and handling space charge

    /** Generate bulk and vacuum meshes using the imported atomistic data */
    int generate_mesh();

    /** Pass mesh to all the objects that need it
     * and transfer temperature from previous iteration to the new mesh */
    int prepare_solvers();

    /** Calculate data of interest on the locations of imported atoms */
    int prepare_export();

    /** Solve Laplace equation on vacuum mesh */
    int solve_laplace(double E0, double V0);

    /** Evolve the PIC simulation one Femocs time step */
    int solve_pic(double advance_time, bool full_run);

    /** Solve transient heat and continuity equations */
    int solve_heat(double T_ambient, double delta_time, bool full_run, int& ccg, int& hcg);

    /** Using the electric field, calculate atomistic charge together with Lorenz and/or Coulomb force */
    int solve_force();

    /** Store the imported atom coordinates and set the flag that enables exporters */
    int finalize(double tstart);

    /** Handle failed subprocess */
    int process_failed(const string &msg);

private:
    /** Determine whether atoms have moved significantly and whether to enable file writing */
    int reinit();

    /** Generate boundary nodes for mesh */
    int generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum);

    /** Transfer mesh from Tetgen into Deal.II */
    int import_mesh();

    /** Interpolate temperature on the centroids of surface quadrangles */
    void calc_surf_temperatures();

    /** Specify mesh address where new mesh will be generated on next run */
    void update_mesh_pointers();

    /** Pick a field solver and calculcate field distribution */
    int run_field_solver();

    /** Pick a heat solver and calculcate temperature & current density distribution */
    int run_heat_solver();

    /** Perform one iteration of PIC calculation */
    int make_pic_step(int& n_lost, int& n_cg, int& n_injected, bool full_run);

    /** Perform one iteration of field emission calculation */
    void calc_heat_emission(bool full_run);

    /** Write restart file so that simulation could be started at t>0 time */
    void write_restart();

    /** Handle Parcas restart file */
    void copy_mdlat();
};

} /* namespace femocs */

#endif /* PROJECTRUNAWAY_H_ */
