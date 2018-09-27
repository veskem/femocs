/*
 * ParticleSpecies.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: andreas, Kyrre, Mihkel
 */

#include "ParticleSpecies.h"

namespace femocs {


ParticleSpecies::ParticleSpecies(double _q_over_m_factor, double _charge, double _Wsp) :
                q_over_m(_q_over_m_factor), q_over_eps0(_charge), Wsp(_Wsp) {}

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

void ParticleSpecies::sort(){
    std::sort(parts.begin(), parts.end());
    //Update the array tracking how many particles per cell
    ordcount.clear();
    if (parts.size() > 0) {
        int cell0 = parts[0].cell;
        ordcount.push_back(0);
        size_t ordCounter = 0;
        for (auto p : parts) {
            if (p.cell == cell0) {
                ordcount[ordCounter] += 1;
            }
            else {
                ordcount.push_back(1);
                cell0 = p.cell;
            }
        }
    }
}

} // namespace femocs
