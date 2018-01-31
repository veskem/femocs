/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas, Mihkel
 */

#include "Pic.h"
#include "Macros.h"

#include <deal.II/base/tensor.h>

namespace femocs {


ParticleSpecies::ParticleSpecies(double _q_over_m_factor, double _charge, double _Wsp) :
                q_over_m_factor(_q_over_m_factor), q_over_eps0(_charge), Wsp(_Wsp) {}

ParticleSpecies::~ParticleSpecies(){
    return;
}

void ParticleSpecies::clear_lost(){
    size_t npart = parts.size();
    size_t nlost = 0;

    //Delete the lost particles from the arrays
    for (size_t i = 0; i < npart; i++) {
        bool islost=false;
        //Is this particle lost?
        for (auto lpart : lost) {
            if (lpart == i) {
                islost=true;
                nlost++;
                break;
            }
        }
        if (nlost==0 or islost) continue; // Don't shuffle this particle left

        parts[i-nlost] = parts[i];
    }

    //Shrink the arrays
    if (nlost > 0){
        parts.resize(npart-nlost);
        cout << "Particles where lost! nlost=" << nlost << endl;
    }

    lost.clear();
}


template<int dim>
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver, fch::CurrentsAndHeating<3> &ch_solver,
        Interpolator &interpolator, EmissionReader &er) : ///< Object to read the temperature data) :
        laplace_solver(laplace_solver), ch_solver(ch_solver), interpolator(interpolator), er(er){}

template<int dim>
Pic<dim>::~Pic() {}

template<int dim>
int Pic<dim>::inject_electrons(const bool fractional_push) {
    vector<Point3> positions, fields;
    vector<int> cells;
    er.inject_electrons(dt, Wsp, positions, fields, cells);

    Point3 velocity0(0,0,0);

    for (int i = 0; i < fields.size(); ++i){
        // Random fractional timestep push -- from a random point [t_(k-1),t_k] to t_k, t_(k + 1/2), using field at t_k.
        if (fractional_push) {
            velocity0 = fields[i] * electrons.q_over_m_factor * dt * ( (double)std::rand()/ RAND_MAX +.5 );
            positions[i] += velocity0 * dt * ( (double)std::rand()/ RAND_MAX );
        }else{
            velocity0 = fields[i] * (electrons.q_over_m_factor * dt * .5);
        }

        // Save to particle arrays
        electrons.inject_particle(positions[i], velocity0, cells[i]);
    }
    return fields.size();

}

//Call the laplace solver with the list of positions and charge(s)
template<int dim>
void Pic<dim>::compute_field(bool first_time) {

    //TODO : save system lhs and copy it . don't rebuild it every time

    double t0;
    start_msg(t0, "=== Solving the Poisson equation...");
    if (first_time)
        laplace_solver.setup_system();

    laplace_solver.assemble_system_lhs();
    laplace_solver.assemble_system_pointcharge(electrons);

    if (anodeBC == "neumann")
        laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top, -E0);
    else if(anodeBC == "dirichlet")
        laplace_solver.assemble_system_dirichlet(fch::BoundaryId::vacuum_top, V0);
    else
        require(false, "ERROR: anodeBC parameter wrong!! anodeBC = " + anodeBC);

    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.);
    laplace_solver.assemble_system_finalize();
    cout << "Solving Poisson. CG iterations = " << laplace_solver.solve() << endl;
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

        //Update the cid_el && check if any particles have left the domain && remove them
        particle.cell = interpolator.update_point_cell(particle.pos, particle.cell);
        if (particle.cell == -1) {
            electrons.lost.push_back(i);
        }
    }
}


template<int dim>
void Pic<dim>::update_fields_and_velocities(){

    //update field
    for (auto particle : electrons.parts) {
        dealii::Point<dim> p(particle.pos.x, particle.pos.y, particle.pos.z);
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(p, particle.cell) ;
        dealii::Point<dim> F(Efield);
        Point3 Field(F);

        //update velocities (corresponds to t + .5dt)
        particle.vel += Field * (dt * electrons.q_over_m_factor) ;

    }
}

template<int dim>
void Pic<dim>::run_cycle(bool first_time) {
    update_positions();
    electrons.sort_parts();
    electrons.clear_lost();
    compute_field(first_time);
    update_fields_and_velocities();
}

//template<int dim>
//void Pic<dim>::clear_lost_particles(){
//    size_t npart = r_el.size();
//    size_t nlost = 0;
//
//    //Delete the lost particles from the arrays
//    for (size_t i = 0; i < npart; i++) {
//        bool islost=false;
//        //Is this particle lost?
//        for (auto lost : lost_el) {
//            if (lost == i) {
//                islost=true;
//                nlost++;
//                break;
//            }
//        }
//        if (nlost==0 or islost) continue; // Don't shuffle this particle left
//
//        r_el[i-nlost] = r_el[i];
//        v_el[i-nlost] = v_el[i];
//        F_el[i-nlost] = F_el[i];
//        cid_el[i-nlost] = cid_el[i];
//    }
//
//    //Shrink the arrays
//    if (nlost > 0){
//        r_el.resize(npart-nlost);
//        v_el.resize(npart-nlost);
//        cid_el.resize(npart-nlost);
//        F_el.resize(npart-nlost);
//        cout << "Particles where lost! nlost=" << nlost << endl;
//    }
//
//    lost_el.clear();
//}

template<int dim>
void Pic<dim>::write_particles(const string filename) {

    // dummy electron always added in the end (to avoid empty electrons crashing ovito)
    ofstream out;
    out.setf(std::ios::scientific);
    out.precision(6);

    string ftype = get_file_type(filename);
    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    cout << "writing particles to " + filename << " n_size = " << electrons.parts.size() << endl;

    out << electrons.parts.size() + 1 << endl;
    out << "Interpolator properties=id:I:1:pos:R:3:vel:R:3:cell:I:1" << endl;

    for (int i = 0; i < electrons.parts.size();  ++i){
        ParticleSpecies::Particle &p = electrons.parts[i];
        out << i << " " << p.pos.x << " " << p.pos.y << " " << p.pos.z << " " <<
        p.vel.x << " " << p.vel.y << " " << p.vel.z << " " << p.cell << endl;
    }

    out << "-1 0.0 0.0 0.0 0.0 0.0 0.0 0" << endl; // always write one dummy particle (ovito crashes with zero parts)

    out.close();
}


//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;
}
