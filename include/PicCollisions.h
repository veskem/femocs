/*
 * PicCollisions.h
 *
 *  Created on: 31.01.2018
 *      Author: Kyrre
 */

#ifndef PICCOLLISIONS_H_
#define PICCOLLISIONS_H_

#include "ParticleSpecies.h"

namespace femocs {

    /*
     * e-e or i-i Coulomb collision routine (for the same type of particles)
     */
    void coll_el_knm_2D( ParticleSpecies &pa);
}

#endif
