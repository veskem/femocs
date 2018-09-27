/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas, Mihkel
 */

#include "Pic.h"
#include "Macros.h"

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>

using namespace std;
namespace femocs {

template<int dim>
Pic<dim>::Pic(PoissonSolver<dim> *poisson, const CurrentHeatSolver<3> *ch_solver,
        const EmissionReader *er, const unsigned int seed) :
        poisson_solver(poisson), ch_solver(ch_solver), emission(er),
        interpolator(er->interpolator),
        electrons(-e_over_me, -e_over_eps0, Wel),
        mersenne{seed}
        {}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions;
    vector<int> cells;
    inject_electrons(positions, cells, *emission->mesh);

    Vec3 velocity(0);

    for (unsigned int i = 0; i < positions.size(); ++i) {
        //update the field
        int cell = interpolator->linhex.deal2femocs(cells[i]);
        Vec3 elfield = interpolator->linhex.interp_gradient(positions[i], cell) ;

        // Random fractional timestep push -- from a random point [t_(k-1),t_k] to t_k, t_(k + 1/2), using field at t_k.
        if (fractional_push) {
            velocity = elfield * (electrons.q_over_m * dt * (uniform(mersenne) + 0.5));
            positions[i] += velocity * (dt * uniform(mersenne));
        } else {
            velocity = elfield * (electrons.q_over_m * dt * 0.5);
        }

        // Save to particle arrays
        electrons.inject_particle(positions[i], velocity, cells[i]);
    }

    inject_stats.injected += positions.size();
    return positions.size();
}

template<int dim>
void Pic<dim>::inject_electrons(vector<Point3> &positions, vector<int> &cells, const TetgenMesh &mesh) {

    const double shift_factor = mesh.tris.stat.edgemin * 1e-6;
    const int n_points = emission->fields->size();

    // loop through quadrangle centroids
    for (int i = 0; i < n_points; ++i) {
        double current = emission->currents[i] * electrons_per_fs;
        double charge = current * dt; //in e
        double n_sps = charge / Wel;

        int intpart = (int) floor(n_sps);
        double frpart = n_sps - intpart;
        int n_electrons = intpart;

        if (uniform(mersenne) < frpart)
            n_electrons++;

        if (n_electrons == 0) continue;

        int quad = abs(emission->fields->get_marker(i));
        int tri = mesh.quads.to_tri(quad);
        int hex = mesh.quad2hex(quad, TYPES.VACUUM);

        hex = interpolator->linhex.femocs2deal(hex);
        Vec3 shift = mesh.tris.get_norm(tri) * shift_factor;

        // generate desired amount of electrons
        // that are uniformly distributed on a given quadrangle
        for (int j = 0; j < n_electrons; ++j) {
            Point3 position = get_rnd_point(quad, mesh);
            // push point little bit inside the vacuum mesh
            position += shift;
            positions.push_back(position);
            cells.push_back(hex);
        }
    }
}

template<int dim>
Point3 Pic<dim>::get_rnd_point(const int quad, const TetgenMesh &mesh) {
    const int tri = mesh.quads.to_tri(quad);
    const int section = quad % n_quads_per_tri;

    int i, j, k;
    if (section == 0) {
        i = 0; j = 1; k = 2;
    } else if (section == 1) {
        i = 1; j = 2; k = 0;
    } else {
        i = 2; j = 0; k = 1;
    }

    SimpleFace sface = mesh.tris[tri];
    Vec3 node0 = mesh.nodes.get_vec(sface[i]);
    Vec3 edge1 = (mesh.nodes.get_vec(sface[j]) - node0) * 0.5;
    Vec3 edge2 = (mesh.nodes.get_vec(sface[k]) - node0) * 0.5;

    array<double,3> bcc;

    // loop until desired point is found
    for (int safe_cntr = 0; safe_cntr < 100; ++safe_cntr) {
        // Generate random point inside parallelogram composed of edge1 & edge2
        double rand1 = uniform(mersenne);
        double rand2 = uniform(mersenne);
        Point3 point = node0 + edge1 * rand1 + edge2 * rand2;

        // calculate barycentric coordinates for a point
        bcc = interpolator->lintri.shape_functions(point, tri);

        // check whether the point is inside the quadrangle
        if (bcc[i] >= bcc[j] && bcc[i] >= bcc[k])
            return point;
    }

    write_silent_msg("Random point generation failed for cell " + d2s(quad));
    return mesh.quads.get_centroid(quad);
}

template<int dim>
int Pic<dim>::update_positions() {

#pragma omp parallel for
    for  (int i = 0; i < electrons.size(); ++i)
        // positions are updated in separate routine to make it easier
        // to parallelize the process with OpenMP
        update_position(i);

    int n_lost_particles = electrons.clear_lost();
    electrons.sort();

    inject_stats.removed += n_lost_particles;
    return n_lost_particles;
}

template<int dim>
void Pic<dim>::update_position(const int particle_index) {
    SuperParticle &particle = electrons.parts[particle_index];

    //update position
    particle.pos += particle.vel * dt;

    //apply periodic boundaries
    // No needed as particles outside the box are deleted anyways
//    particle.pos.x = periodic_image(particle.pos.x, box.xmax, box.xmin);
//    particle.pos.y = periodic_image(particle.pos.y, box.ymax, box.ymin);

    // Update the cell ID; if any particles have left the domain their ID is set to -1
    // and they will be removed once we call clear_lost
    const bool b1 = particle.pos.x > box.xmin && particle.pos.x < box.xmax;
    const bool b2 = particle.pos.y > box.ymin && particle.pos.y < box.ymax;
    const bool b3 = particle.pos.z < box.zmax;

    if (b1 && b2 && b3)
        particle.cell = update_point_cell(particle);
    else
        particle.cell = -1;
}

template<int dim>
int Pic<dim>::update_point_cell(const SuperParticle& particle) const {
    int femocs_cell = interpolator->linhex.deal2femocs(particle.cell);
    femocs_cell = interpolator->linhex.locate_cell(particle.pos, femocs_cell);
    if (femocs_cell < 0) return -1;
    return interpolator->linhex.femocs2deal(femocs_cell);
}

template<int dim>
void Pic<dim>::update_velocities(){
    for (auto &particle : electrons.parts) {
        // find electric field
        int cell = interpolator->linhex.deal2femocs(particle.cell);
        Vec3 elfield = interpolator->linhex.interp_gradient(particle.pos, cell);

        // update velocities (corresponds to t + .5dt)
        particle.vel += elfield * (dt * electrons.q_over_m);
    }
}

template<int dim>
void Pic<dim>::collide_particles() {
    if (!coll_coulomb_ee) return;

    // TODO masses and charges
    double variance_factor = dt * lanlog / (twopi*eps0*eps0);

    vector<vector<size_t>> particles_in_cell;
    group_and_shuffle_particles(particles_in_cell);

    for (size_t cell = 0; cell < particles_in_cell.size(); ++cell) {
        size_t n_parts = particles_in_cell[cell].size();

        // constant factor in Coulomb collisions for this cell
        // TODO why is volume needed ?
        double Acoll_cell = variance_factor * n_parts / poisson_solver->get_cell_vol(cell);

        // TODO what happens with unpaired electrons?
        for (size_t p = 0; p < n_parts; p+=2) {
            size_t p1 = particles_in_cell[cell][p];
            size_t p2 = particles_in_cell[cell][p+1];

            // relative velocity
            Vec3 v_rel = electrons.parts[p1].vel - electrons.parts[p2].vel;
            double v_rel_norm  = v_rel.norm();

            // if u->0, the relative change might be big,
            // but it's a big change on top of nothing => SKIP
            if (v_rel_norm < 1.e-20) continue;

            // MC to find the scattering angle, through transformation of theta
            // which is is distributed with P(theta) = Gauss(0, <theta^2>)
            double variance2 = Acoll_cell / (v_rel_norm*v_rel_norm*v_rel_norm);

            // TODO sqrt or not?
            std::normal_distribution<double> rnd_gauss(0.0, sqrt(variance2));
            double theta = rnd_gauss(mersenne);     // scattering angle
            double phi = twopi * uniform(mersenne); // azimuth angle

            // calculate rotation/scattering matrix entries
            double sin_theta = sin(theta);
            double cos_phi = cos(phi);
            double cos_theta_minus_one = cos(theta) - 1;

            double v_cross  = sqrt(v_rel.x * v_rel.x + v_rel.y * v_rel.y);
            double A = v_rel_norm * sin_theta * sin(phi) / v_cross;
            double B = v_rel.x * sin_theta * cos_phi / v_cross;
            double C = v_rel.y * sin_theta * cos_phi / v_cross;

            // calculate change of velocity
            Vec3 v_delta = Vec3(
                    v_rel.dotProduct( Vec3(cos_theta_minus_one, A, B) ),
                    v_rel.dotProduct( Vec3(-A, cos_theta_minus_one, C) ),
                    v_rel.dotProduct( Vec3(-B, -C, cos_theta_minus_one) )
            );
            v_delta *= 0.5;

            // Update the particle velocities (identical particle masses)
            electrons.parts[p1].vel += v_delta;
            electrons.parts[p2].vel -= v_delta;
        } // loop over particles
    } // loop over cells
}

template<int dim>
void Pic<dim>::group_and_shuffle_particles(vector<vector<size_t>> &parts_in_cell) {
    const int n_cells = interpolator->linhex.size();
    const int n_particles = electrons.parts.size();

    // group particles
    parts_in_cell = vector<vector<size_t>>(n_cells);
    for (size_t p = 0; p < n_particles; ++p) {
        int cell = interpolator->linhex.deal2femocs(electrons.parts[p].cell);
        if (cell >= 0) parts_in_cell[cell].push_back(p);
    }

    // randomize the ordering of particles in a group
    for (int cell = 0; cell < n_cells; ++cell) {
        int n_particles = parts_in_cell[cell].size();
        if (n_particles > 1) {
            std::uniform_int_distribution<int> rnd_int(0, n_particles-1);

            for (size_t &p : parts_in_cell[cell]) {
                int new_index = rnd_int(mersenne);
                size_t temporary = parts_in_cell[cell][new_index];
                parts_in_cell[cell][new_index] = p;
                p = temporary;
            }
        }
    }
}

template<int dim>
void Pic<dim>::collide_particles_old() {
    if (!coll_coulomb_ee) return;

//    static const double Amplcoulomb =  1.;  // Amplification of coulomb collisions, only for testing purposes
//    static const double LanLog = 13.;       // Coulomb Log
//    Constant factor in coulomb collisions
//    double Acoll = Amplcoulomb * ( !kind ? 1.0 : dt_ion*dt_ion*dt_ion ) * LanLog * SQU(SQU(Omega_pe)) * ncoll /
//      ( TWOPI*PI * SQU(SQU(dz))*SQU(dz) * SQU(pa->mass) * SQU(Ndb) * N_sp );

    vector<size_t> inds2coll; //Not nice for parallelization (neither over spieces- or cell)
    double Acoll = 1.0; //TODO !!!!!

    size_t Next = 0; // Index of the first particle in the current cell
    // (the way this is done today messes up paralellization)
    for (auto ord : electrons.ordcount) {
        //Number of particles left to collide in the current cell
        size_t N2coll = ord;
        //Resize the indices
        inds2coll.resize(N2coll);

        //Constant factor in columb collisions for this cell
        double Acoll_cell = Acoll * ord / poisson_solver->get_cell_vol(electrons.parts[Next].cell);

        if ( N2coll>1 ) {
            //Pick two random particle indices (j,k) from the particles that have not yet collide
            for (size_t loopind  = 0; loopind < N2coll; loopind++ ) {
                inds2coll[loopind]=loopind;
            }
            size_t j(0), k(0);
            while ( N2coll>0 ) {
                if ( N2coll>1 ) {
                    size_t jind = (size_t)(N2coll * uniform(mersenne));
                    if (jind == N2coll) jind--;
                    j = inds2coll[jind];
                    for (size_t loopind = jind + 1; loopind < N2coll; loopind++) {
                        inds2coll[loopind-1] = inds2coll[loopind];
                    }
                    N2coll--;

                    jind = (size_t)(N2coll * uniform(mersenne));
                    if (jind == N2coll) jind--;
                    k = inds2coll[jind];
                    for (size_t loopind = jind + 1; loopind < N2coll; loopind++) {
                        inds2coll[loopind-1] = inds2coll[loopind];
                    }
                    N2coll--;
                }
                else {
                    //Last particle in a cell with odd number of particles:
                    // collide it with "j" from last time
                    k = inds2coll[N2coll-1];
                    N2coll--;
                }

                // Collide the pair as described in:
                //   K. Matyash,
                //   Kinetic Modeling of Multi-Component Edge Plasmas,
                //   PhD thesis, Ernst-Moritz-Arndt-Universität Greifswald, 2003
                //   http://home.mpcdf.mpg.de/~knm/thesis/2.pdf
                // and
                //   Takizuka and H. Abe,
                //   A binary collision model for plasma simulation with a particle code,
                //   Journal of Computational Physics 25 (1977) 205

                // relative velocity
                Vec3 v_rel = electrons.parts[j+Next].vel - electrons.parts[k+Next].vel;
                // W = ||v_rel||
                double W  = v_rel.norm();
                if (W < 1.e-10) {
                    //if u->0, the relative change might be big,
                    // but it's a big change on top of nothing => SKIP
                    continue;
                }

                //MC to find the scattering angle, through transformation of delta
                // which is is distributed with P(delta) = N(0, <delta^2>)
                double deltaVar2 = Acoll_cell/(W*W*W); // Variance <delta^2>
                //double delta = GausRandom(0,sqrt(deltaVar2)); // TODO !!!!!!
                double delta = (uniform(mersenne) - 0.5) * deltaVar2;  // TODO !!!!!! Temporary hack until we have GausRandom
                double sinTheta = 2 * delta / (1+delta*delta);
                double oneMinusCosTheta = sinTheta*delta;
                //Scattering azimuth angle
                double Phi = twopi * uniform(mersenne);
                double cosPhi = cos(Phi);
                double sinPhi = sin(Phi);

                //v_rel projected into XY plane (and inverse)
                Vec3 v_rel_DELTA;
                double v_rel_XY  = sqrt(v_rel.x * v_rel.x + v_rel.y * v_rel.y);
                if (v_rel_XY > 1e-10) {
                    //Normal case, need rotation of coordinate system etc.
                    v_rel_DELTA.x =   (v_rel.x/v_rel_XY)*v_rel.z * sinTheta*cosPhi
                                    - (v_rel.y/v_rel_XY)*W * sinTheta*sinPhi
                                    - v_rel.x * oneMinusCosTheta;
                    v_rel_DELTA.y =   (v_rel.y/v_rel_XY)*v_rel.z * sinTheta*cosPhi
                                    + (v_rel.x/v_rel_XY)*W * sinTheta*sinPhi
                                    -  v_rel.y * oneMinusCosTheta;
                    v_rel_DELTA.z = - v_rel_XY*sinTheta*cosPhi
                                    - v_rel.z*oneMinusCosTheta;
                }
                else {
                    //v_rel in t-direction only
                    // => coordinate systems are identical
                    // => simpler equations
                    v_rel_DELTA.x = W*sinTheta*cosPhi;
                    v_rel_DELTA.y = W*sinTheta*sinPhi;
                    v_rel_DELTA.z = -W*oneMinusCosTheta;
                }

                //Update the particle velocities (identical particle masses)
                electrons.parts[j + Next].vel += v_rel_DELTA/2.;
                electrons.parts[k + Next].vel -= v_rel_DELTA/2.;
            }//END do-loop over particles
        }//END if (N2coll > 1)
        Next += ord;
    }//END loop over cells
}

template<int dim>
void Pic<dim>::write(const string &filename) const {
    if (!MODES.WRITEFILE) return;

    string ftype = get_file_type(filename);
    require(ftype == "xyz" || ftype == "movie", "Invalid file type: " + ftype);

    ofstream out;
    out.setf(std::ios::scientific);
    out.precision(6);

    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    out << max(1, electrons.size()) << endl;
    out << "time= " << GLOBALS.TIME << ", Pic properties=id:I:1:pos:R:3:Velocity:R:3:cell:I:1" << endl;

    out.setf(ios::scientific);
    out.precision(6);

	// Ovito can't handle 0 particles, so in case of empty system write dummy particle
    if (electrons.size() == 0) {
		out << "-1 0.0 0.0 0.0 0.0 0.0 0.0 0" << endl;
	} else {
		for (int i = 0; i < electrons.size();  ++i)
		    out << i << " " << electrons.parts[i] << endl;
	}

    out.close();
}

//Tell the compiler which types to actually compile, so that they are available for the linker
template class Pic<3>;

} // namespace femocs
