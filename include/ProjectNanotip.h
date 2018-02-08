/*
 * ProjectNanotip.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef PROJECTNANOTIP_H_
#define PROJECTNANOTIP_H_

using namespace std;
namespace femocs {

#include "AtomReader.h"
#include "Config.h"
#include "SolutionReader.h"
#include "Surface.h"
#include "Interpolator.h"
#include "TetgenMesh.h"
#include "laplace.h"
#include "GeneralProject.h"

/*
 *
 */
class ProjectNanotip : public GeneralProject {
public:

    ProjectNanotip(const AtomReader &reader, const Config &conf);

    ~ProjectNanotip() {}

    int run();

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const double elfield, const int timestep=-1);

    /** Function to export the calculated electric field on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param Ex        x-component of electric field
     * @param Ey        y-component of electric field
     * @param Ez        z-component of electric field
     * @param Enorm     norm of electric field
     * @return          success of the operation (always 0)
     */
    int export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Function to linearly interpolate electric field at points anywhere in space
     * @param n_points  number of points where electric field is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Function to linearly interpolate electric field at points on or near the surface.
     * It is faster than interpolate_elfield but is inadequate for points far away from surface
     * @param n_points  number of points where electric field is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Function to linearly interpolate electric potential at given points
     * @param n_points  number of points where electric potential is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param phi       electric potential
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_phi(const int n_points, const double* x, const double* y, const double* z, double* phi, int* flag);

    /** Generate bulk and vacuum meshes using the imported atomistic data */
    int generate_mesh();

    /** Solve Laplace equation on vacuum mesh */
    int solve_laplace(const double E0);

    /** Determine whether atoms have moved significantly and whether to enable file writing */
    int reinit(const int timestep);

    /** Store the imported atom coordinates and set the flag that enables exporters */
    int finalize();

    /** Force the data to the files for debugging purposes */
    int force_output();

protected:
    static constexpr double delta_t_MD = 4.05e-15; ///< MD timestep in seconds

    bool skip_meshing;      ///< If the mesh is to be kept the same
    bool fail;              ///< If some process failed
    double t0;              ///< CPU timer
    int timestep;           ///< counter to measure how many times Femocs has been called
    int last_full_timestep; ///< last time step Femocs did full calculation
    string timestep_string; ///< time step written to file name
    vector<int> atom2face;  ///< surface atom to triangle index map

    Surface dense_surf;     ///< non-coarsened surface atoms
    Surface extended_surf;  ///< atoms added for the surface atoms

    Interpolator vacuum_interpolator = Interpolator("elfield", "potential");
    fch::Laplace<3> laplace_solver;   ///< Laplace equation solver

    /** Generate boundary nodes for mesh */
    int generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum);
};

} /* namespace femocs */

#endif /* PROJECTNANOTIP_H_ */
