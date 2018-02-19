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

class ParticleSpecies {
public:
    ParticleSpecies(double q_ovr_m, double charge, double Wsp);
    ~ParticleSpecies() {};

    void inject_particle(const Point3 &pos, const Vec3 &vel, const int cell) {
        parts.push_back(SuperParticle(pos, vel, cell));
    }

    int clear_lost();

    void sort();

    int size() const { return parts.size(); }

    const double q_over_m_factor; ///< charge/mass [A^2 / (V fs^2)]
    const double q_over_eps0;     ///< (whole) particle charge / eps0 [e/VÅ]
    const double Wsp;             ///< SP weight [particles/superparticle]

    vector<SuperParticle> parts;
    vector<size_t> ordcount;
};

} // namespace femocs

#endif /* PARTICLESPECIES_H_ */
