/*
 * ProjectRunaway.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef PROJECTRUNAWAY_H_
#define PROJECTRUNAWAY_H_

#include "GeneralProject.h"
#include "ProjectNanotip.h"
#include "AtomReader.h"
#include "Config.h"
#include "SolutionReader.h"
#include "Surface.h"
#include "Interpolator.h"
#include "TetgenMesh.h"
#include "CurrentsAndHeating.h"
#include "CurrentsAndHeatingStationary.h"
#include "PhysicalQuantities.h"
#include "Pic.h"

using namespace std;
namespace femocs {

/**
 * Class for calculating electric field and heating effects
 * around a nanostructure that fits into central cylinder
 * by solving Laplace or Poisson equation
 * and taking into account change of temperature due to field emission.
 */
class ProjectRunaway : public ProjectNanotip {
public:

    ProjectRunaway(AtomReader &a, Config &c);
    ~ProjectRunaway() {}

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const double elfield, const int timestep=-1);

    /** Evolve the PIC simulation one Femocs time step */
    int solve_pic(const double E0, const double dt_main);

    /** Pick a method to solve heat and continuity equations on bulk mesh */
    int solve_heat(const double T_ambient);

    /** Force the data to the files for debugging purposes */
    int force_output();

private:
    Interpolator bulk_interpolator = Interpolator("rho", "temperature");

    FieldReader surface_fields;       ///< fields on surface hex face centroids
    HeatReader  surface_temperatures; ///< temperatures & current densities on surface hex face centroids
    EmissionReader emission;          ///< emission data on centroids of surface quadrangles

    /// physical quantities used in heat calculations
    fch::PhysicalQuantities phys_quantities = fch::PhysicalQuantities(conf.heating);
    fch::CurrentsAndHeatingStationary<3>  ch_solver1;     ///< first    steady-state currents and heating solver
    fch::CurrentsAndHeatingStationary<3>  ch_solver2;     ///< second   steady-state currents and heating solver
    fch::CurrentsAndHeatingStationary<3>* ch_solver;      ///< active   steady-state currents and heating solver
    fch::CurrentsAndHeatingStationary<3>* prev_ch_solver; ///< previous steady-state currents and heating solver
    fch::CurrentsAndHeating<3> ch_transient_solver;       ///< transient currents and heating solver

    Pic<3> pic_solver;    ///< Class for solving Poisson equation and handling space charge

    /** Solve steady-state heat and continuity equations */
    int solve_stationary_heat();

    /** Solve transient heat and continuity equations */
    int solve_transient_heat(const double delta_time);

    /** Solve transient heat and continuity equation until convergence is reached */
    int solve_converge_heat();

    /** Import meshes to dealii and set params to various objects */
    int prepare_fem();
};

} /* namespace femocs */

#endif /* PROJECTRUNAWAY_H_ */
