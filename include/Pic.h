/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas, Mihkel
 */

#ifndef PIC_H_
#define PIC_H_

#include "Interpolator.h"
#include "EmissionReader.h"
#include "Config.h"
#include "TetgenMesh.h"
#include "TetgenCells.h"
#include "Primitives.h"
#include "PoissonSolver.h"
#include "CurrentHeatSolver.h"
#include "ParticleSpecies.h"

#include <random>

namespace femocs {

/** Class for running PIC (particle-in-cell) simulations to solve Poisson equation.
 * For further details about PIC, see Kyrre Ness Sjøbæk PhD thesis at
 * https://cds.cern.ch/record/2226840
 */
template<int dim>
class Pic {
public:
    Pic(PoissonSolver<dim> *poisson_solver, const CurrentHeatSolver<3> *ch_solver,
            const EmissionReader *emission, const unsigned int seed);
    ~Pic() {};

    /** Inject electrons according to the field emission surface distribution */
    int inject_electrons(const bool fractional_push);

    /** e-e or i-i Coulomb collision routine (for the same type of particles) */
    void collide_particles();
    void collide_particles_old();

    /** Update the positions of the particles and the cell they belong. */
    int update_positions();

    /** Update the velocities on the particles using the fields calculated at the given positions */
    void update_velocities();

    /** Write the particle data in the current state in movie file */
    void write(const string &filename) const;
    
    /** Store various data */
    void set_params(const Config::PIC &conf, const double _dt, const TetgenNodes::Stat _box) {
        dt = _dt;
        Wel = conf.Wsp_el;
        electrons.set_Wsp(Wel);
        box = _box;
        coll_coulomb_ee = conf.coll_coulomb_ee;
        lanlog = conf.lanlog;
    }

    /** Return pointer to charged super-particles */
    ParticleSpecies* get_particles() { return &electrons; }

    int get_n_electrons() const { return electrons.size(); }

    void stats_reinit(){
        inject_stats.injected = 0;
        inject_stats.removed = 0;
    }

    void reinit(){
        stats_reinit();
        electrons.clear();
    }

    int get_injected(){ return inject_stats.injected; }

    int get_removed(){ return inject_stats.removed; }

    /**  Check if the number of injected particles is roughly equal to the number of removed*/
    bool is_stable(){
        return abs(inject_stats.injected - inject_stats.removed) / inject_stats.injected < 0.1 ;
    }

private:
    /// charge/mass for electrons for multiplying the velocity update
    static constexpr double e_over_me = 17.58820024182468;
    static constexpr double e_over_eps0 = 180.9512268; ///< electron charge / vacuum permittivity [V*Angstrom]
    static constexpr double electrons_per_fs = 6.2415e3; ///< definition of 1 ampere
    static constexpr double twopi = 6.2831853071795864;  ///< 2 * pi
    static constexpr double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]

    //Parameters
    double Wel = .01;   ///< Super particle weighting [particles/superparticle]
    double dt = 1.;     ///< timestep [fs]
    bool coll_coulomb_ee = false;   ///< Switch 2e->2e Coulomb collisions on/off
    double lanlog = 0;  ///< Landau logarithm

    PoissonSolver<dim> *poisson_solver;        ///< object to solve Poisson equation in the vacuum mesh
    const CurrentHeatSolver<dim> *ch_solver;  ///< transient currents and heating solver
    const EmissionReader *emission;         ///< object to obtain the field emission data
    const Interpolator *interpolator;       ///< data & operation for interpolating in vacuum
    TetgenNodes::Stat box;                  ///< simubox size data

    ParticleSpecies electrons;   ///< super particles carring charge density in space

    struct InjectStats {
        int injected = 0;
        int removed = 0;
    } inject_stats;

    mt19937 mersenne;     ///< Mersenne twister pseudo-random number engine
    uniform_real_distribution<double> uniform{0.0, 1.0};  ///< Random nr generator that maps Mersenne twister output uniformly into range [0.0 1.0] (both inclusive)

    /** Clear all particles out of the box */
    void clear_lost_particles();

    /** Computes the charge density for each FEM DOF */
    void compute_field(bool first_time = false, bool write_time = false);

    /** Inject electron super particles at the surface faces, depending on the current and the timestep */
    void inject_electrons(vector<Point3> &pos, vector<int> &cells, const TetgenMesh &mesh);

    /** Generate point with an uniform distribution inside a quadrangle */
    Point3 get_rnd_point(const int quad, const TetgenMesh &mesh);

    /** Update position and cell index of a super particle with given index */
    void update_position(const int particle_index);

    /** Find the hexahedron where the point is located.
     * Both input and output hex indices are Deal.II ones. */
    int update_point_cell(const SuperParticle& particle) const;

    void group_and_shuffle_particles(vector<vector<size_t>> &particles_in_cell);
};

} // namespace femocs

#endif /* PIC_H_ */
