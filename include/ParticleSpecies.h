/*
 * ParticleSpecies.h
 *
 *  Created on: Jan 30, 2018
 *      Author: andreas, Mihkel
 */

#ifndef PARTICLESPECIES_H_
#define PARTICLESPECIES_H_

#include "Primitives.h"
#include <algorithm>

namespace femocs {

/** Super particles for PIC simulation */
class ParticleSpecies {
public:
    ParticleSpecies(double q_over_m, double q_over_eps0, double Wsp);
    ~ParticleSpecies() {};

    void inject_particle(const Point3 &pos, const Vec3 &vel, const int cell) {
        parts.push_back(SuperParticle(pos, vel, cell));
    }

    void clear() { parts.clear(); }

    int clear_lost();

    int size() const { return parts.size(); }

    double get_Wsp() const { return Wsp; }
    void set_Wsp(double Wsp) {
        require(Wsp >= 0, "Invalid SP weight: " + d2s(Wsp));
        this->Wsp = Wsp;
    }

    /** Accessor for accessing the i-th SP */
    const SuperParticle& operator [](const size_t i) const {
        require(i < size(), "Invalid index: " + d2s(i));
        return parts[i];
    }
    SuperParticle& operator [](const size_t i) {
        require(i < size(), "Invalid index: " + d2s(i));
        return parts[i];
    }

    /** Iterator to access the particles */
    typedef Iterator<ParticleSpecies, SuperParticle> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    const double q_over_m;      ///< SP charge / its mass [A^2 / (V fs^2)]
    const double q_over_eps0;   ///< (whole) particle charge / eps0 [e/VÅ]

private:
    double Wsp;                 ///< SP weight [particles/superparticle]
    vector<SuperParticle> parts;
};

} // namespace femocs

#endif /* PARTICLESPECIES_H_ */
