/*
 * ParticleSpecies.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: andreas, Kyrre, Mihkel
 */

#include "ParticleSpecies.h"

namespace femocs {

ParticleSpecies::ParticleSpecies(double _q_over_m, double _q_over_eps0, double _Wsp) :
        q_over_m(_q_over_m), q_over_eps0(_q_over_eps0), Wsp(_Wsp)
{}

int ParticleSpecies::clear_lost() {
    size_t npart = parts.size();
    size_t nlost = 0;

    // Delete the lost particles from the array
    for (size_t i = 0; i < npart; i++) {
        if (parts[i].cell==-1)
            nlost++;
        else if (nlost > 0)
            parts[i-nlost] = parts[i];
    }

    // Shrink the array
    if (nlost > 0) parts.resize(npart-nlost);
    return nlost;
}

} // namespace femocs
