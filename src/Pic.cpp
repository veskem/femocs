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
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver, fch::CurrentsAndHeating<3> &ch_solver,
        Interpolator &interpolator, EmissionReader &er) :
        laplace_solver(laplace_solver), ch_solver(ch_solver), interpolator(interpolator), er(er),
        coll_coulomb_ee(false) {}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions;
    vector<Vec3> fields;
    vector<int> cells;
    er.inject_electrons(dt, Wsp, positions, fields, cells);

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
void Pic<dim>::compute_field(bool first_time) {

    //TODO : save system lhs and copy it . don't rebuild it every time

    double t0;
    start_msg(t0, "\nSetting up system lhs... first_time = " + to_string(first_time));
    laplace_solver.setup_system(first_time);
    if (first_time){
        laplace_solver.assemble_system_lhs();
        laplace_solver.save_system();
    }else{
        laplace_solver.restore_system();
    }
    end_msg(t0);

    start_msg(t0, "\nSetting up system RHS and BCs...");
    laplace_solver.assemble_system_pointcharge(electrons);
    if (anodeBC == "neumann")
        laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top, -E0);
    else if(anodeBC == "dirichlet")
        laplace_solver.assemble_system_dirichlet(fch::BoundaryId::vacuum_top, V0);
    else
        require(false, "ERROR: anodeBC parameter wrong!! anodeBC = " + anodeBC);

    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.);
    laplace_solver.assemble_system_finalize();
    end_msg(t0);

    start_msg(t0, "\nSolving Poisson equation... CG iterations = ");
    cout << laplace_solver.solve() << " ";
    end_msg(t0);
}

template<int dim>
void Pic<dim>::update_positions(){

    for  (size_t i = 0; i < electrons.parts.size(); ++i){
        ParticleSpecies::Particle &particle = electrons.parts[i];

        //update position
        particle.pos += particle.vel * dt;

        //apply periodic boundaries
        particle.pos.x = periodic_image(particle.pos.x, box.xmax, box.xmin);
        particle.pos.y = periodic_image(particle.pos.y, box.ymax, box.ymin);

        //Update the cellID; if any particles have left the domain their ID is set to -1
        // and they will be removed once we call clear_lost
        particle.cell = interpolator.update_point_cell(particle.pos, particle.cell);
    }
}

template<int dim>
void Pic<dim>::update_velocities(){
    //update field
    for (auto particle : electrons.parts) {
        dealii::Point<dim> p(particle.pos.x, particle.pos.y, particle.pos.z);
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(p, particle.cell) ;

        //update velocities (corresponds to t + .5dt)
        particle.vel += Vec3(Efield) * (dt * electrons.q_over_m_factor) ;
    }
}

template<int dim>
void Pic<dim>::run_cycle(bool first_time) {
    double t0;

    start_msg(t0,"=== Updating PIC particle positions ...");
    update_positions();
    electrons.clear_lost();
    electrons.sort_parts();
    end_msg(t0);

    start_msg(t0, "=== Calculating electric field ...");
    compute_field(first_time);
    end_msg(t0);

    start_msg(t0, "=== Updating the PIC particle velocities ...");
    update_velocities();
    end_msg(t0);

    start_msg(t0, "=== Colliding PIC particles ...");
    if (coll_coulomb_ee){
        coll_el_knm_2D(electrons, dt, laplace_solver);
    }
    end_msg(t0);
}

template<int dim>
void Pic<dim>::write_particles(const string filename, double time) {
    string ftype = get_file_type(filename);
    require(ftype == "xyz" || ftype == "movie", "Invalid file type: " + ftype);

    ofstream out;
    out.setf(ios::scientific);
    out.precision(6);

    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    out << electrons.parts.size() + 1 << endl;
    out << "time = " << time << " Interpolator properties=id:I:1:pos:R:3:vel:R:3:cell:I:1" << endl;

    for (int i = 0; i < electrons.parts.size();  ++i){
        ParticleSpecies::Particle &p = electrons.parts[i];
        out << i << " " << p.pos << " " << p.vel << " " << p.cell << endl;
    }

    out << "-1 0.0 0.0 0.0 0.0 0.0 0.0 0" << endl; // always write one dummy particle (ovito crashes with zero parts)

    out.close();
}

//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;

} // namespace femocs
