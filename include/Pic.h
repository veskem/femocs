/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#ifndef PIC_H_
#define PIC_H_

#include "laplace.h"
#include "mesh_preparer.h"
#include "SolutionReader.h"
#include "currents_and_heating.h"

#include <deal.II/base/point.h>

namespace femocs {


template<int dim>
class Pic {
public:
    Pic(fch::Laplace<dim> &laplace_solver, FieldReader &fr,
            fch::CurrentsAndHeating<3> &ch_solver, HeatReader &hr, EmissionReader &er);
    ~Pic();

    /**Injects electrons
     *     Indexing: (x1 y1 [z1] x2 y2 [z2] ...)
     */
    int injectElectrons(const double* const r, const size_t n);

    /**
     * Runs a full pic cycle : inject, push, update field
     */

    /**
     * Inject electrons according to the field emission surface distribution
     */
    int injectElectrons();

    /**
     * Run an particle and field update cycle
     */
    void runCycle();

    /**
     * Write the particle data in the current state in movie file
     */
    void writeParticles(const string filename);
    
    void set_timestep(double _dt){
        dt = _dt;
    }

    void set_E0(double _E0){
        E0 = _E0;
    }




private:

    //ELECTRONS

    std::vector<dealii::Point<dim>> r_el; ///< Particle positions [Å]

    std::vector<dealii::Point<dim>> v_el; ///< Particle velocities [Å/fs = 10 um/s]

    std::vector<dealii::Point<dim>> F_el; ///< Particle forces [eV / Å]
    //Management
    std::vector<int> cid_el; //Index of the cell where the particle is inside
    std::vector<int> lost_el; //Index into the cell array containing lost particles
    
    //std::vector<double> charges; // charges

    //Constants
    const double q_over_m_factor = 17.58820024182468; // 1[e]/511e3[eV/c**2]*E[V/Å]*dt[fs] = 5866.792...*E*dt [Å/fs] charge/mass for electrons for multiplying the velocity update
    const double q_over_eps0 = 180.9512268; // particle charge [e] / epsilon_0 [e/VÅ] = 1 [e] * (8.85...e-12/1.6...e-19/1e10 [e/VÅ])**-1
    
    const double Wsp = .01; // Super particle weighting (particles/superparticle)

    double dt = 1.; //timestep

    double E0 = -1; //Applied field at Neumann boundary


    fch::Laplace<dim> &laplace_solver; ///< Laplace solver object to solve the Poisson in the vacuum mesh
    fch::CurrentsAndHeating<3> &ch_solver;       ///< transient currents and heating solver
    FieldReader &fr; ///< Object to read the electric field
    HeatReader &hr; ///< Object to read the temperature data
    EmissionReader &er; ///< Object to calculate the emission data

    /**
     * Clear all particles out of the box
     */
    void clearLostParticles();

    /**
     * Update the positions of the particles and the cell they belong.
     */

    void updatePositions();

    void updateFieldAndVelocities();


    /**Computes the charge density for each FEM DOF
     *
     */
    void computeField();
};

}

#endif
