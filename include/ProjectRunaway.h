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

    int run(const int timestep=-1);

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const double elfield, const int timestep=-1);

    /** Force the data to the files for debugging purposes */
    int force_output();

    /** Export the solution on the imported atomistic points.
     * @param n_points  # of first imported atoms where data exported; 0 disables the export
     * @param cmd   label of the data to be exported. See Labels class for a list of possible cmd-s
     * @param data  array where results are written. Vector data is exported coordinate-wise, i.e in a form x1,y1,z1,x2,y2,...
     */
    int export_results(const int n_points, const string &cmd, double* data);

    /** Export the solution on the specified points.
      * @param n_points    # of first imported atoms where data exported; 0 disables the export
      * @param cmd         label of the data to be exported. See Labels class for a list of possible cmd-s
      * @param on_surface  data is located on or near the surface.
      * @param x,y,z       coordinates of the points where interpolation is performed
      * @param data        array where results are written. Vector data is exported coordinate-wise, i.e in a form x1,y1,z1,x2,y2,...
      * @param flag        array showing whether specified point was located inside the mesh (1) or not (0)
      */
    int interpolate_results(const int n_points, const string &cmd, const bool on_surface,
            const double* x, const double* y, const double* z, double* data, int* flag);

private:
    bool fail;                  ///< If some process failed
    double t0;                  ///< CPU timer
    int timestep;               ///< counter to measure how many times Femocs has been called
    double last_write_time = -1.e200; ///< Keeps the time that last file output was done
    double last_heat_time = -1.e200; ///< Last time heat was updated

    bool mesh_changed = false;          ///< True if new mesh has been created
    int last_full_timestep; ///< last time step Femocs did full calculation
    string timestep_string; ///< time step written to file name
    // as surface atom->triangle mapping is quite heavy but useful in many places, it's good to prevent doing it many times
    vector<int> atom2face;  ///< surface atom to triangle index map

    Interpolator vacuum_interpolator;  ///< data & operations for interpolating field & potential in vacuum
    Interpolator bulk_interpolator;    ///< data & operations for interpolating current density & temperature in bulk

    Surface dense_surf;       ///< non-coarsened surface atoms
    Surface extended_surf;    ///< atoms added for the surface atoms

    FieldReader surface_fields;       ///< fields on surface hex face centroids
    HeatReader  surface_temperatures; ///< temperatures & current densities on surface hex face centroids

    PhysicalQuantities phys_quantities; ///< quantities used in heat calculations
    PoissonSolver<3> poisson_solver;    ///< Poisson equation solver
    CurrentHeatSolver<3> ch_solver;     ///< transient currents and heating solver

    EmissionReader emission;          ///< emission data on centroids of surface quadrangles
    Pic<3> pic_solver;                       ///< class for solving Poisson equation and handling space charge

    /** Write output data to files */
    int write_results();

    /** Check if enough time has passed since the last file write_results */
    bool write_time() const { return GLOBALS.TIME < (last_write_time + conf.behaviour.write_period); }

    /** Determine whether atoms have moved significantly and whether to enable file writing */
    int reinit(const int timestep);

    /** Store the imported atom coordinates and set the flag that enables exporters */
    int finalize();

    /** Generate boundary nodes for mesh */
    int generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum);

    /** Generate bulk and vacuum meshes using the imported atomistic data */
    int generate_mesh();

    /** Read bulk and vacuum meshes from file and generate metadata for them */
    int read_mesh();

    /** Solve Laplace equation on vacuum mesh */
    int solve_laplace(double E0);

    /** Evolve the PIC simulation one Femocs time step */
    int solve_pic(const double advance_time);

    int converge_pic(double max_time);

    /** Solve transient heat and continuity equations */
    int solve_heat(double T_ambient, double delta_time, int& ccg, int& hcg);

    /** Solve transient heat and continuity equation until convergence is reached */
    int converge_heat(double T_ambient);
    
    /** Import mesh to FEM solvers and initialize interpolators */
    int prepare_solvers();

    /** Calculate data of interest on the locations of imported atoms */
    int prepare_export();
};

} /* namespace femocs */

#endif /* PROJECTRUNAWAY_H_ */
