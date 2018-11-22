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
Pic<dim>::Pic(const PoissonSolver<dim> *poisson, const EmissionReader *er,
        const Interpolator *i, const Config::PIC *config, const unsigned int seed) :
        poisson_solver(poisson), emission(er), interpolator(i), conf(config),
        electrons(-e_over_me, -e_over_eps0, 0),
        mersenne{seed}
{}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions;
    vector<int> cells;
    inject_electrons(positions, cells, *emission->mesh);

    if (positions.size() > conf->max_injected)
        return -1 * positions.size();

    Vec3 velocity(0);

    for (unsigned int i = 0; i < positions.size(); ++i) {
        //update the field
        int cell = interpolator->linhex.deal2femocs(cells[i]);
        Vec3 elfield = interpolator->linhex.interp_gradient(positions[i], cell) ;

        // Random fractional timestep push -- from a random point [t_(k-1),t_k] to t_k, t_(k + 1/2), using field at t_k.
        if (fractional_push) {
            velocity = elfield * (electrons.q_over_m * data.dt * (uniform(mersenne) + 0.5));
            positions[i] += velocity * (data.dt * uniform(mersenne));
        } else {
            velocity = elfield * (electrons.q_over_m * data.dt * 0.5);
        }

        // Save to particle arrays
        electrons.inject_particle(positions[i], velocity, cells[i]);
    }

    data.injected += positions.size();
    return positions.size();
}

template<int dim>
void Pic<dim>::inject_electrons(vector<Point3> &positions, vector<int> &cells, const TetgenMesh &mesh) {

    const double shift_factor = mesh.tris.stat.edgemin * 1e-6;
    const int n_points = emission->fields->size();

    // loop through quadrangle centroids
    for (int i = 0; i < n_points; ++i) {
        double current = emission->currents[i] * electrons_per_fs;
        double charge = current * data.dt; //in e
        double n_sps = charge / electrons.get_Wsp();

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

    data.removed += n_lost_particles;
    return n_lost_particles;
}

template<int dim>
void Pic<dim>::update_position(const int particle_index) {
    SuperParticle &electron = electrons[particle_index];

    //update position
    electron.pos += electron.vel * data.dt;

    bool b1 = true;
    bool b2 = true;
    bool b3 = electron.pos.z < data.box.zmax;

    if (conf->periodic) {
        // apply periodic boundaries
        electron.pos.x = periodic_image(electron.pos.x, data.box.xmax, data.box.xmin);
        electron.pos.y = periodic_image(electron.pos.y, data.box.ymax, data.box.ymin);
    } else {
        // check the boundaries in x,y-direction
        b1 = electron.pos.x > data.box.xmin && electron.pos.x < data.box.xmax;
        b2 = electron.pos.y > data.box.ymin && electron.pos.y < data.box.ymax;
        if (!b1 || !b2) {
            write_silent_msg("Electron " + d2s(particle_index) + " crossed "
                    "simubox x or y boundary and will be deleted.\n  "
                    "Consider increasing simubox width or making PIC periodic.");
        }
    }

    // Update the cell ID; if any particles have left the domain their ID is set to -1
    // and they will be removed once we call clear_lost
    if (b3 && b2 && b1)
        electron.cell = update_point_cell(electron);
    else
        electron.cell = -1;
}

template<int dim>
int Pic<dim>::update_point_cell(const SuperParticle& particle) const {
    // in case mesh has changed, the particle.cell has quite random value w.r the new mesh.
    // However the cell nr from previous mesh is little bit better guess than just 0,
    // as there is a hope, that new mesh was generated similarly to the old one
    // and therefore also the cell indices in old and new mesh are similar although different.
    int femocs_cell = interpolator->linhex.deal2femocs(particle.cell);
    femocs_cell = interpolator->linhex.locate_cell(particle.pos, femocs_cell);
    if (femocs_cell < 0) return -1;
    return interpolator->linhex.femocs2deal(femocs_cell);
}

template<int dim>
void Pic<dim>::update_velocities(){
    for (SuperParticle &electron : electrons) {

        // find electric field
        int cell = interpolator->linhex.deal2femocs(electron.cell);
        Vec3 elfield = interpolator->linhex.interp_gradient(electron.pos, cell);

        // update velocities (corresponds to t + .5dt)
        electron.vel += elfield * (data.dt * electrons.q_over_m);
    }
}

/* Collide same kind of charged particles as described in
 * Takizuka and H. Abe
 * A binary collision model for plasma simulation with a particle code
 * Journal of Computational Physics 25 (1977) 205 */
template<int dim>
void Pic<dim>::collide_particles() {
    if (!conf->collide_ee) return;

    double variance_factor = e_over_me * e_over_eps0 * electrons.get_Wsp();
    variance_factor *= variance_factor * data.dt * conf->landau_log / twopi;

    vector<vector<size_t>> particles_in_cell;
    group_and_shuffle_particles(particles_in_cell);

    int cell = 0;
    for (vector<size_t> &particles : particles_in_cell) {
        size_t n_parts = particles.size();
        if (n_parts <= 1) continue;

        // constant factor in Coulomb collisions for this cell
        double Acoll_cell = variance_factor * n_parts / poisson_solver->get_cell_vol(cell++);
        size_t start_p = 0;

        // for odd number of particles, perform first three collisions differently
        if (n_parts % 2) {
            start_p = 3;
            size_t p1 = particles[0];
            size_t p2 = particles[1];
            size_t p3 = particles[2];

            collide_pair(p1, p2, 0.5*Acoll_cell);
            collide_pair(p1, p3, 0.5*Acoll_cell);
            collide_pair(p2, p3, 0.5*Acoll_cell);
        }

        // perform pair-wise collision
        for (size_t p = start_p; p < n_parts; p+=2) {
            size_t p1 = particles[p];
            size_t p2 = particles[p+1];
            collide_pair(p1, p2, Acoll_cell);
        }
    }
}

template<int dim>
void Pic<dim>::collide_pair(int p1, int p2, double variance_factor) {
    // Relative velocity
    Vec3 v_rel = electrons[p1].vel - electrons[p2].vel;
    double v_rel_norm  = v_rel.norm();

    // If u->0, the relative change might be big,
    // but it's a big change on top of nothing => SKIP
    if (v_rel_norm < 1.e-10) return;

    // Use MC to find the scattering angle, through transformation of delta
    // which is is distributed with P(delta) = Gauss(0, <delta^2>)
    double variance2 = variance_factor / (v_rel_norm*v_rel_norm*v_rel_norm);

    std::normal_distribution<double> rnd_gauss(0.0, sqrt(variance2));
    double delta = rnd_gauss(mersenne);     // scattering angle
    double phi = twopi * uniform(mersenne); // azimuth angle

    // Calculate rotation/scattering matrix entries
    double sin_theta = 2 * delta / (1 + delta*delta);
    double cos_theta_minus_one = -delta * sin_theta;
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double v_perp  = sqrt(v_rel.x * v_rel.x + v_rel.y * v_rel.y);

    Vec3 v_delta;
    if (v_perp > 1e-10) {
        /* Here we prefer to use the matrix in a form given by
         * Tskhakaya et al
         * The Particle‐In‐Cell Method
         * Contributions to Plasma Physics, 47(8‐9), pp.563-594, 2007
         *
         * It is computationally equivalent to the one by Takizuka and Abe,
         * but reveals better the symmetry.
         */
        double a = v_rel_norm * sin_theta * sin_phi / v_perp;
        double bx = v_rel.x * sin_theta * cos_phi / v_perp;
        double by = v_rel.y * sin_theta * cos_phi / v_perp;

        v_delta.x = v_rel.dotProduct( Vec3(cos_theta_minus_one, a, bx) );
        v_delta.y = v_rel.dotProduct( Vec3(-a, cos_theta_minus_one, by) );
        v_delta.z = v_rel.dotProduct( Vec3(-bx, -by, cos_theta_minus_one) );
    } else {
        v_delta.x = v_rel_norm * sin_theta * cos_phi;
        v_delta.y = v_rel_norm * sin_theta * sin_phi;
        v_delta.z = v_rel_norm * cos_theta_minus_one;
    }
    v_delta *= 0.5;

    // Update the particle velocities (particle masses are identical)
    electrons[p1].vel += v_delta;
    electrons[p2].vel -= v_delta;
}

template<int dim>
void Pic<dim>::group_and_shuffle_particles(vector<vector<size_t>> &parts_in_cell) {
    const int n_cells = poisson_solver->get_n_cells();
    const int n_particles = electrons.size();

    // Group particles
    parts_in_cell = vector<vector<size_t>>(n_cells);
    for (size_t p = 0; p < n_particles; ++p) {
        int cell = electrons[p].cell;
        require(cell < n_cells, "Invalid cell " + d2s(cell) + " associated with superparticle " + d2s(p));
        if (cell >= 0) parts_in_cell[cell].push_back(p);
    }

    // Randomize the ordering of particles in a group
    for (vector<size_t> &particles : parts_in_cell) {
        int n_particles = particles.size();
        // no need to shuffle <= 3 particles, as their pairing doesn't change
        if (n_particles > 3) {
            std::uniform_int_distribution<int> rnd_int(0, n_particles-1);

            for (size_t &p : particles) {
                int new_index = rnd_int(mersenne);
                std::swap(p, particles[new_index]);
            }
        }
    }
}

template<int dim>
void Pic<dim>::read(const string &filename) {
    string ftype = get_file_type(filename);
    require(ftype == "bin", "Unimplemented file type: " + ftype);

    ifstream in(filename);
    require(in.is_open(), "Can't open a file " + file_name);

    string str;
    while (in >> str) {
        if (str == "$Electrons") {
            int n_electrons;
            in >> n_electrons >> GLOBALS.TIME >> GLOBALS.TIMESTEP;
            getline(in, str);

            electrons.reserve(n_electrons);

            SuperParticle electron;
            for (int i = 0; i < n_electrons; ++i) {
                in.read(reinterpret_cast<char*>(&electron), sizeof(SuperParticle));
                electrons.inject_particle(electron);
            }
        }
    }

    in.close();
}

template<int dim>
void Pic<dim>::write(const string &filename) const {
    if (!MODES.WRITEFILE) return;

    string ftype = get_file_type(filename);
    ofstream out;
    if (ftype == "movie" || ftype == "bin")
        out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    if (ftype == "xyz" || ftype == "movie")
        write_xyz(out);
    else if (ftype == "bin")
        write_bin(out);
    else
        require(false, "Invalid file type: " + ftype);

    out.close();
}

template<int dim>
void Pic<dim>::write_xyz(ofstream &out) const {
    out << max(1, electrons.size()) << endl;
    out << "Time=" << GLOBALS.TIME << ", Timestep=" << GLOBALS.TIMESTEP
            << ", Pic properties=id:I:1:pos:R:3:Velocity:R:3:cell:I:1" << endl;

    out.setf(ios::scientific);
    out.precision(6);

    // Ovito can't handle 0 particles, so in case of empty system write dummy particle
    if (electrons.size() == 0) {
        out << "-1 0.0 0.0 0.0 0.0 0.0 0.0 0" << endl;
    } else {
        for (int i = 0; i < electrons.size();  ++i)
            out << i << " " << electrons[i] << endl;
    }
}

template<int dim>
void Pic<dim>::write_bin(ofstream &out) const {
    const int n_data = electrons.size();

    // write data header
    out << "$Electrons\n"
            << n_data << " " << GLOBALS.TIME << " " << GLOBALS.TIMESTEP << "\n";

    // write data
    for (SuperParticle const &sp : electrons) {
        out.write ((char*)&sp, sizeof (SuperParticle));
    }

    // close data header
    out << "\n$EndElectrons\n";
}

//Tell the compiler which types to actually compile, so that they are available for the linker
template class Pic<3>;

} // namespace femocs
