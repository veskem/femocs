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
#include "FileWriter.h"

#include <random>

namespace femocs {

/** Class for running PIC (particle-in-cell) simulations to solve Poisson equation.
 * For further details about PIC, see Kyrre Ness Sjøbæk PhD thesis at
 * https://cds.cern.ch/record/2226840
 */
template<int dim>
class Pic: public FileWriter {
public:
    Pic(const PoissonSolver<dim> *poisson_solver, const EmissionReader *emission,
            const Interpolator *interpolator, const Config::PIC *conf, const unsigned int seed);
    ~Pic() {};

    /** Inject electrons according to the field emission surface distribution */
    int inject_electrons(const bool fractional_push);

    /** e-e or i-i Coulomb collision routine (for the same type of particles) */
    void collide_particles();

    /** Update the positions of the particles and the cell they belong. */
    int update_positions();

    /** Update the velocities on the particles using the fields calculated at the given positions */
    void update_velocities();

    /** Read particle data from file */
    void read(const string &filename);
    
    /** Store various data */
    void set_params(const double dt, const TetgenNodes::Stat box) {
        electrons.set_Wsp(conf->weight_el);
        data.dt = dt;
        data.box = box;
    }

    /** Return pointer to charged super-particles */
    ParticleSpecies* get_particles() { return &electrons; }

    /** Return number of stored electron super particles */
    int get_n_electrons() const { return electrons.size(); }

    /** Return entry for xyz file */
    int size() const { return max(1, electrons.size()); }

    void stats_reinit() {
        data.injected = 0;
        data.removed = 0;
    }

    void reinit() {
        stats_reinit();
        electrons.clear();
    }

    /** Return # electron super particles that were injected during last push */
    int get_injected() const { return data.injected; }

    /** Return # electron super particles that were deleted during last clean */
    int get_removed() const { return data.removed; }

    /**  Check if the # injected particles is roughly equal to the # removed */
    bool is_stable() const {
        return abs(data.injected - data.removed) / data.injected < 0.1 ;
    }

private:
    static constexpr double e_over_me = 17.58820024182468;///< charge / mass of electrons for multiplying the velocity update
    static constexpr double e_over_eps0 = 180.9512268;    ///< electron charge / vacuum permittivity [V*Angstrom]
    static constexpr double electrons_per_fs = 6.2415e3;  ///< definition of 1 ampere
    static constexpr double twopi = 6.2831853071795864;   ///< 2 * pi

    const PoissonSolver<dim> *poisson_solver; ///< object to solve Poisson equation in the vacuum mesh
    const EmissionReader *emission;           ///< object to obtain the field emission data
    const Interpolator *interpolator;         ///< data & operation for interpolating in vacuum
    const Config::PIC *conf;                  ///< PIC configuration parameters

    ParticleSpecies electrons;    ///< Electron super particles carrying charge density in space

    struct Data {
        TetgenNodes::Stat box;    ///< simubox size data
        double dt = 0;            ///< timestep [fs]
        int injected = 0;
        int removed = 0;
    } data;

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

    /** Group SP-s by the common cells and shuffle their ordering before performing Coulomb collision */
    void group_and_shuffle_particles(vector<vector<size_t>> &particles_in_cell);

    /** Perform binary collision between two charged particles of the same kind */
    void collide_pair(int p1, int p2, double variance);

    /** Specify file types that can be written */
    bool valid_extension(const string &ext) const {
        return ext == "xyz" || ext == "movie" || ext == "restart";
    }

    /** Label for restart data */
    string get_restart_label() const { return "Electrons"; }

    /** Write the particle data into xyz file */
    void write_xyz(ofstream &out) const;

    /** Write the particle data into restart file */
    void write_bin(ofstream &out) const;
};

} // namespace femocs

#endif /* PIC_H_ */
