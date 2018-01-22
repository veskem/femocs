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
    Pic(/*fch::Laplace<dim> &laplace_solver, fch::CurrentsAndHeating<3> &ch_solver,
            FieldReader &fr, HeatReader &hr, EmissionReader &er*/);
    ~Pic();

    /**Injects electrons
     *     Indexing: (x1 y1 [z1] x2 y2 [z2] ...)
     */
    int injectElectrons(const double* const r, const size_t n, FieldReader &fr);

    int injectElectrons();

    /**Computes the charge density for each FEM DOF
     *
     */
    void computeField(const double E0);

    /** Pushes the particles given the fields for a delta time dt [sec]
     *
     */
    void pushParticles(const double dt, FieldReader &fr);

    /** Write the position and velocities of the particles to a file
     *
     */
    void writeParticles(const string filename);
    
private:

    //ELECTRONS

    std::vector<dealii::Point<dim>> r_el; //Particle positions [Å]
    std::vector<dealii::Point<dim>> v_el; //Particle velocities [Å/fs = 10 um/s]
    std::vector<int> cid_el; ///< Index of the cell where the particle is inside

    std::vector<double> charges; ///< charges

    //Constants
    const double q_over_m = 1.0; ///< [?] charge/mass for electrons
    const double q_over_eps0 = 180.9512268; ///< particle charge [e] / epsilon_0 [e/VA]
    const double W = 1.0; ///< Super particle weighting


    fch::Laplace<dim> &laplace_solver; ///< Laplace solver object to solve the Poisson in the vacuum mesh
    fch::CurrentsAndHeating<3> &ch_solver;       ///< transient currents and heating solver
    FieldReader &fr; ///< Object to read the electric field
    HeatReader &hr; ///< Object to read the temperature data
    EmissionReader &er; ///< Object to calculate the emission data
};

}

#endif
