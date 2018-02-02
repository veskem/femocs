/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#ifndef PIC_H_
#define PIC_H_

#include "laplace.h"
//#include "mesh_preparer.h"
#include "Interpolator.h"
#include "currents_and_heating.h"
#include "SolutionReader.h"
#include "Config.h"
#include "TetgenCells.h"
#include "Primitives.h"
#include "ParticleSpecies.h"
#include "PicCollisions.h"


#include <deal.II/base/point.h>
#include <algorithm>

namespace femocs {

template<int dim>
class Pic {
public:
    Pic(fch::Laplace<dim> &laplace_solver, fch::CurrentsAndHeating<3> &ch_solver, Interpolator &interpolator, EmissionReader &er);
    ~Pic();

    /**
     * Inject electrons according to the field emission surface distribution
     */
    int inject_electrons(const bool fractional_push);

    /**
     * Run an particle and field update cycle
     */
    void run_cycle(bool first_time = false);

    /**
     * Write the particle data in the current state in movie file
     */
    void write_particles(const string filename, double time);
    
    void set_params(const Config::Field &conf_lap,
            const Config::PIC &conf_pic,
            double _dt, TetgenNodes::Stat _box){
        dt = _dt;
        Wsp = conf_pic.Wsp_el;
        E0 = conf_lap.E0;
        V0 = conf_lap.V0;
        anodeBC = conf_lap.anodeBC;
        box = _box;
        coll_coulomb_ee = conf_pic.coll_coulomb_ee;
    }

private:


    //Constants
    const double e_over_m_e_factor = 17.58820024182468; ///< charge/mass for electrons for multiplying the velocity update
    const double e_over_eps0 = 180.9512268; ///< particle charge [e] / epsilon_0 [e/VÅ] = 1 [e] * (8.85...e-12/1.6...e-19/1e10 [e/VÅ])**-1
    
    //Parameters
    double Wsp = .01;   ///< Super particle weighting [particles/superparticle]
    double dt = 1.;     ///< timestep [fs]
    double E0 = -1;     ///< Applied field at Neumann boundary [V/Å]
    double V0 = 1000;   ///< Applied voltage (in case Dirichlet boundary at anode) [V]
    string anodeBC = "Neumann"; ///< Boundary type at the anode


    fch::Laplace<dim> &laplace_solver;      ///< Laplace solver object to solve the Poisson in the vacuum mesh
    fch::CurrentsAndHeating<3> &ch_solver;  ///< transient currents and heating solver
    Interpolator &interpolator;
    EmissionReader &er;                     ///< Object to calculate the emission data

    TetgenNodes::Stat box;                  ///< Object containing box data

    bool coll_coulomb_ee;                   ///< Switch 2e->2e Coulomb collisions on/off
    
    ParticleSpecies electrons = ParticleSpecies(-e_over_m_e_factor, -e_over_eps0, Wsp);


    /** Clear all particles out of the box */
    void clear_lost_particles();

    /** Update the positions of the particles and the cell they belong. */
    void update_positions();

    /** Update the velocities on the particles using the fields calculated at the given positions */
    void update_velocities();

    /** Computes the charge density for each FEM DOF */
    void compute_field(bool first_time = false);
    
};

}

#endif
