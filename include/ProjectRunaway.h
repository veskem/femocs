/*
 * ProjectRunaway.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef PROJECTRUNAWAY_H_
#define PROJECTRUNAWAY_H_

using namespace std;
namespace femocs {

#include "GeneralProject.h"
#include "ProjectNanotip.h"
#include "AtomReader.h"
#include "Config.h"
#include "SolutionReader.h"
#include "Surface.h"
#include "Interpolator.h"
#include "TetgenMesh.h"
#include "physical_quantities.h"
#include "currents_and_heating.h"
#include "currents_and_heating_stationary.h"
#include "Pic.h"

/*
 *
 */
class ProjectRunaway : public ProjectNanotip {

    ProjectRunaway(const AtomReader &a, const Config &c);
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

    /** Function to export the calculated temperatures on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param T         temperature in the atom location
     * @return          success of the operation (always 0)
     */
    int export_temperature(const int n_atoms, double* T);

    /** Calculate and export charges & forces on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param xq        charges and forces in PARCAS format (xq[0] = q1, xq[1] = Fx1, xq[2] = Fy1, xq[3] = Fz1, xq[4] = q2, xq[5] = Fx2 etc)
     * @return          success of the operation (always 0)
     */
    int export_charge_and_force(const int n_atoms, double* xq);

    /** Export Laplace + Coulomb force and pair potential on imported atoms
     * @param n_atoms  number of first atoms the data will be exported
     * @param xnp      forces in PARCAS format & units (xnp[0] = Fx1, xnp[1] = Fy1, xnp[2] = Fz1, xnp[3] = Fx2 etc)
     * @param Epair    potential energy per atom
     * @param Vpair    total potential energy of atoms. Pot. due to Coloumb forces are added here. NOTE: Lorentz is missing!
     * @return         success of the operation (always 0)
     */
    int export_force_and_pairpot(const int n_atoms, double* xnp, double* Epair, double* Vpair);

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
