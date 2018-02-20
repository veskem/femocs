/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas, Mihkel
 */

#include "Pic.h"
#include "Macros.h"

#include <deal.II/base/tensor.h>

using namespace std;
namespace femocs {

template<int dim>
Pic<dim>::Pic(fch::PoissonSolver<dim> *poisson, fch::CurrentHeatSolver<3> *ch_solver,
        Interpolator *interpolator, EmissionReader *er) :
        poisson_solver(poisson), ch_solver(ch_solver), interpolator(interpolator), er(er),
        coll_coulomb_ee(false) {}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions;
    vector<Vec3> fields;
    vector<int> cells;
    er->inject_electrons(dt, Wsp, positions, fields, cells);

    Vec3 velocity(0);

    for (int i = 0; i < fields.size(); ++i){
        // Random fractional timestep push -- from a random point [t_(k-1),t_k] to t_k, t_(k + 1/2), using field at t_k.
        if (fractional_push) {
            velocity = fields[i] * electrons.q_over_m_factor * dt * ( (double)rand()/ RAND_MAX +.5 );
            positions[i] += velocity * dt * ( (double)rand()/ RAND_MAX );
        } else {
            velocity = fields[i] * (electrons.q_over_m_factor * dt * .5);
        }

        // Save to particle arrays
        electrons.inject_particle(positions[i], velocity, cells[i]);
    }
    return fields.size();
}

//Call the laplace solver with the list of positions and charge(s)
template<int dim>
int Pic<dim>::compute_field(bool first_time) {
    poisson_solver->assemble_poisson(first_time);
    return poisson_solver->solve();
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
    //update field
    for (auto particle : electrons.parts) {
        dealii::Point<dim> p(particle.pos.x, particle.pos.y, particle.pos.z);
        dealii::Tensor<1,dim> Efield = poisson_solver->probe_efield(p, particle.cell) ;

        //update velocities (corresponds to t + .5dt)
        particle.vel += Vec3(Efield) * (dt * electrons.q_over_m_factor) ;
    }
}

template<int dim>
int Pic<dim>::run_cycle(bool first_time) {
    // calculate electric field
    int n_cg_steps = compute_field(first_time);

    // updating the PIC particle velocities
    update_velocities();

    // collide PIC particles
    if (coll_coulomb_ee)
        coll_el_knm_2D(electrons, dt, *poisson_solver);

    return n_cg_steps;
}

template<int dim>
void Pic<dim>::write(const string filename, const double time) const {
    string ftype = get_file_type(filename);
    require(ftype == "xyz" || ftype == "movie", "Invalid file type: " + ftype);

    if (electrons.size() == 0) return;  // Ovito can't handle 0 particles

    ofstream out;

    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    out << electrons.size() << endl;
    out << "time=" << time << ", Pic properties=id:I:1:pos:R:3:Velocity:R:3:cell:I:1" << endl;

    out.setf(ios::scientific);
    out.precision(6);

    for (int i = 0; i < electrons.size();  ++i)
        out << i << " " << electrons.parts[i] << endl;

    out.close();
}

//Tell the compiler which types to actually compile, so that they are available for the linker
template class Pic<3>;

} // namespace femocs
