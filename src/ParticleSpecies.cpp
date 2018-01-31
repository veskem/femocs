/*
 * ParticleSpecies.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: andreas, Kyrre
 */

#include "ParticleSpecies.h"

namespace femocs {


ParticleSpecies::ParticleSpecies(double _q_over_m_factor, double _charge, double _Wsp) :
                q_over_m_factor(_q_over_m_factor), q_over_eps0(_charge), Wsp(_Wsp) {}

ParticleSpecies::~ParticleSpecies(){
    return;
}

void ParticleSpecies::clear_lost(){
    size_t npart = parts.size();
    size_t nlost = 0;

    //Delete the lost particles from the arrays
    for (size_t i = 0; i < npart; i++) {
        //If this particle is lost or nothing has been lost:
        // skip shuffling this particle
        if (parts[i].cell==-1) {
            nlost++;
            continue;
        }
        else if (nlost==0) {
            continue;
        }

        parts[i-nlost] = parts[i];
    }
    
    //Shrink the arrays
    if (nlost > 0){
        parts.resize(npart-nlost);
        cout << "Particles were lost! nlost=" << nlost << endl;
    }

}

void ParticleSpecies::sort_parts(){
    sort(parts.begin(), parts.end());
    //Update the array tracking how many particles per cell
    ordcount.clear();
    if (parts.size() > 0) {
        int cell0 = parts[0].cell;
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

}
