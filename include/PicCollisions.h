/*
 * PicCollisions.h
 *
 *  Created on: 31.01.2018
 *      Author: Kyrre
 */

#ifndef PICCOLLISIONS_H_
#define PICCOLLISIONS_H_

#include "ParticleSpecies.h"
#include "PoissonSolver.h" //Needs to access the mesh

namespace femocs {

    /*
     * e-e or i-i Coulomb collision routine (for the same type of particles)
     */
    void coll_el_knm_2D( ParticleSpecies &pa, const double dt, PoissonSolver<3> &poisson_solver );
}

#endif
