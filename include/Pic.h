/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas, Mihkel
 */

#ifndef PIC_H_
#define PIC_H_

#include "Interpolator.h"
#include "SolutionReader.h"
#include "Config.h"
#include "TetgenCells.h"
#include "Primitives.h"
#include "ParticleSpecies.h"
#include "PicCollisions.h"
#include "DealSolver.h"
#include <deal.II/base/point.h>
#include <algorithm>
#include "CurrentHeatSolver.h"
#include "PoissonSolver.h"

namespace femocs {

template<int dim>
class Pic {
public:
    Pic(fch::PoissonSolver<dim> *poisson_solver, fch::CurrentHeatSolver<3> *ch_solver,
            Interpolator *interpolator, EmissionReader *er);
    ~Pic() {};

    /** Inject electrons according to the field emission surface distribution */
    int inject_electrons(const bool fractional_push);

    /**
     * Run an particle and field update cycle
     */
    int run_cycle(bool first_time = false, bool write_time = false);

    /** Update the positions of the particles and the cell they belong. */
    int update_positions();

    /** Write the particle data in the current state in movie file */
    void write(const string filename, const double time) const;
    
    /** Store various data */
    void set_params(const Config::Field &conf_lap, const Config::PIC &conf_pic, const double _dt,
            const TetgenNodes::Stat _box) {
        dt = _dt;
        Wel = conf_pic.Wsp_el;
        electrons.set_Wsp(Wel);
        E0 = conf_lap.E0;
        V0 = conf_lap.V0;
        anodeBC = conf_lap.anodeBC;
        box = _box;
        coll_coulomb_ee = conf_pic.coll_coulomb_ee;
    }

    /** Return pointer to charged super-particles */
    ParticleSpecies* get_particles() { return &electrons; }

private:
    /// charge/mass for electrons for multiplying the velocity update
    static constexpr double e_over_m_e_factor = 17.58820024182468;
    /// particle charge [e] / epsilon_0 [e/VÅ] = 1 [e] * (8.85e-12 / 1.6e-19/1e10 [e/VÅ])**-1
    static constexpr double e_over_eps0 = 180.9512268;

    //Parameters
    double Wel = .01;   ///< Super particle weighting [particles/superparticle]
    double dt = 1.;     ///< timestep [fs]
    double E0 = -1;     ///< Applied field at Neumann boundary [V/Å]
    double V0 = 1000;   ///< Applied voltage (in case Dirichlet boundary at anode) [V]
    string anodeBC = "Neumann"; ///< Boundary type at the anode
    bool coll_coulomb_ee;       ///< Switch 2e->2e Coulomb collisions on/off

    fch::PoissonSolver<dim> *poisson_solver; ///< object to solve Poisson equation in the vacuum mesh
    fch::CurrentHeatSolver<3> *ch_solver;    ///< transient currents and heating solver
    Interpolator *interpolator;
    EmissionReader *er;                     ///< object to calculate the emission data
    TetgenNodes::Stat box;                  ///< simubox size data
    
    ParticleSpecies electrons = ParticleSpecies(-e_over_m_e_factor, -e_over_eps0, Wel);


    /** Clear all particles out of the box */
    void clear_lost_particles();

    /** Update the velocities on the particles using the fields calculated at the given positions */
    void update_velocities();

    /** Computes the charge density for each FEM DOF */
    void compute_field(bool first_time = false, bool write_time = false);
    
};

} // namespace femocs

#endif /* PIC_H_ */
