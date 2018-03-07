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

/** Class for running PIC (particle-in-cell) simulations.
 * For further details about PIC, see Kyrre Ness Sjøbæk PhD thesis at
 * https://cds.cern.ch/record/2226840
 */
template<int dim>
Pic<dim>::Pic(fch::PoissonSolver<dim> *poisson, const fch::CurrentHeatSolver<3> *ch_solver,
        const EmissionReader *er, const unsigned int seed) :
        poisson_solver(poisson), ch_solver(ch_solver), emission(er), interpolator(er->interpolator),
        mersenne{seed}
        {}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions;
    vector<int> cells;
    inject_electrons(positions, cells, *emission->mesh);

    Vec3 velocity(0);

    for (int i = 0; i < positions.size(); ++i) {
        //update the field
        dealii::Point<dim> p(positions[i].x, positions[i].y, positions[i].z);
        dealii::Tensor<1,dim> Efield = poisson_solver->probe_efield(p, cells[i]) ;

        // Random fractional timestep push -- from a random point [t_(k-1),t_k] to t_k, t_(k + 1/2), using field at t_k.
        if (fractional_push) {
            velocity = Vec3(Efield) * (electrons.q_over_m_factor * dt * (uniform(mersenne) + 0.5));
            positions[i] += velocity * (dt * uniform(mersenne));
        } else {
            velocity = Vec3(Efield) * (electrons.q_over_m_factor * dt * 0.5);
        }

        // Save to particle arrays
        electrons.inject_particle(positions[i], velocity, cells[i]);
    }

    return positions.size();
}

template<int dim>
int Pic<dim>::inject_electrons(vector<Point3> &positions, vector<int> &cells,
        const TetgenMesh &mesh) {

    const double shift_factor = mesh.tris.stat.edgemin * 1e-5;
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
    const int section = quad % mesh.quads.n_quads_per_tri;

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
        interpolator->lintri.get_shape_functions(bcc, point, tri);

        // check whether the point is inside the quadrangle
        if (bcc[i] >= bcc[j] && bcc[i] >= bcc[k])
            return point;
    }

    write_silent_msg("Random point generation failed for cell " + to_string(quad));
    return mesh.quads.get_centroid(quad);
}

template<int dim>
int Pic<dim>::update_positions() {
    for  (size_t i = 0; i < electrons.size(); ++i) {
        SuperParticle &particle = electrons.parts[i];

        //update position
        particle.pos += particle.vel * dt;

        //apply periodic boundaries
        particle.pos.x = periodic_image(particle.pos.x, box.xmax, box.xmin);
        particle.pos.y = periodic_image(particle.pos.y, box.ymax, box.ymin);

        // Update the cell ID; if any particles have left the domain their ID is set to -1
        // and they will be removed once we call clear_lost
        particle.cell = interpolator->update_point_cell(particle.pos, particle.cell);
    }

    int n_lost_particles = electrons.clear_lost();
    electrons.sort();
    return n_lost_particles;
}

template<int dim>
void Pic<dim>::update_velocities(){

    for (auto &particle : electrons.parts) {
        //get field
        dealii::Point<dim> p(particle.pos.x, particle.pos.y, particle.pos.z);
        dealii::Tensor<1,dim> Efield = poisson_solver->probe_efield(p, particle.cell) ;

        //update velocities (corresponds to t + .5dt)
        particle.vel += Vec3(Efield) * (dt * electrons.q_over_m_factor) ;
    }
}

template<int dim>
int Pic<dim>::run_cycle(bool first_time, bool charge_write_time) {
    // calculate electric field
    poisson_solver->assemble_poisson(first_time, charge_write_time);
    int n_cg_steps = poisson_solver->solve();

    // updating the PIC particle velocities
    update_velocities();

    // collide PIC particles
    if (coll_coulomb_ee)
        coll_el_knm_2D(electrons, dt, *poisson_solver);

    return n_cg_steps;
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
    out << "time=" << GLOBALS.TIME << ", Pic properties=id:I:1:pos:R:3:Velocity:R:3:cell:I:1" << endl;

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
