/*
 * ProjectNanotip.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef PROJECTNANOTIP_H_
#define PROJECTNANOTIP_H_

#include "GeneralProject.h"
#include "AtomReader.h"
#include "Config.h"
#include "SolutionReader.h"
#include "Surface.h"
#include "Interpolator.h"
#include "laplace.h"

using namespace std;
namespace femocs {

/**
 * Class for calculating electric field around a nanostructure that fits into central cylinder
 * by solving Laplace equation.
 */
class ProjectNanotip : public GeneralProject {
public:

    ProjectNanotip(AtomReader &reader, Config &conf);
    ~ProjectNanotip() {};

    int run(const int timestep=-1);

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const double elfield, const int timestep=-1);

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

    int export_results(const int n_points, const string &cmd, double* data);

    int interpolate_results(const int n_points, const string &cmd, const bool surface,
            const double* x, const double* y, const double* z, double* data, int* flag);

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
