/*
 * ParticleSpecies.h
 *
 *  Created on: Jan 30, 2018
 *      Author: andreas, Mihkel
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

    void inject_particle(Point3 &_pos, Vec3 &_vel, int _cell){
        parts.push_back(Particle{.pos = _pos, .vel = _vel, .cell = _cell});
    }

    void clear_lost();

    struct Particle{
        Point3 pos; ///< Particle position [Å]
        Vec3 vel; ///< Particle velocity [Å/fs]
        int cell;   ///< Cell ID, -1 if outside the domain
                    ///< (to be cleared up by ParticleSpecies::clear_lost())
        bool operator < (const Particle &partj) {return (cell < partj.cell);}
    };

    void sort_parts();

    const double q_over_m_factor; ///< charge/mass [A^2 / (V fs^2)]
    const double q_over_eps0;     ///< (whole) particle charge / eps0 [e/VÅ]
    const double Wsp;             ///< SP weight [particles/superparticle]

    vector<Particle> parts;
    vector<size_t> ordcount;
};

} // namespace femocs

#endif /* PARTICLESPECIES_H_ */
