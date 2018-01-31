/*
 * ParticleSpecies.h
 *
 *  Created on: Jan 30, 2018
 *      Author: andreas
 */

#ifndef PARTICLESPECIES_H_
#define PARTICLESPECIES_H_

#include "Primitives.h"

#include <deal.II/base/point.h>
#include <algorithm>

namespace femocs {


class ParticleSpecies{
public:
    ParticleSpecies(double q_ovr_m, double charge, double Wsp);
    ~ParticleSpecies();

    void inject_particle(Point3 _pos, Point3 _vel, int _cell){
        parts.push_back(Particle{.pos = _pos, .vel = _vel, .cell = _cell});
    }

    void clear_lost();

    struct Particle{
        Point3 pos;
        Point3 vel;
        int cell;
        bool operator < (const Particle &partj) {return (cell < partj.cell);}
    };

    void sort_parts();

    const double q_over_m_factor; ///< charge/mass [A^2 / (V fs^2)]
    const double q_over_eps0; ///< charge / eps0
    const double Wsp; ///< SP weight

    vector<Particle> parts;
    vector<size_t> ordcount;
};

} // namespace




#endif /* INCLUDE_PARTICLESPECIES_H_ */
