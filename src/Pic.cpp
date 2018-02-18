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

/* ==================================================================
 *  ========================== PicSolver ===========================
 * ================================================================== */
 
template<int dim>
PicSolver<dim>::PicSolver(fch::PoissonSolver<dim> &poisson_solver, Interpolator &interpolator, EmissionReader &er) :
        poisson_solver(poisson_solver), interpolator(interpolator), er(er),
        coll_coulomb_ee(false) {}

template<int dim>
int PicSolver<dim>::inject_electrons(const bool fractional_push) {
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
int PicSolver<dim>::compute_field(bool first_time) {

    //TODO : save system lhs and copy it . don't rebuild it every time

    // Set-up system lhs
    poisson_solver.setup_system(first_time);
    if (first_time){
        poisson_solver.assemble_system_lhs();
        poisson_solver.save_system();
    } else
        poisson_solver.restore_system();

    // Set-up system RHS and BCs
    poisson_solver.assemble_system_pointcharge(electrons);
    if (anodeBC == "neumann")
        poisson_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top, -E0);
    else if(anodeBC == "dirichlet")
        poisson_solver.assemble_system_dirichlet(fch::BoundaryId::vacuum_top, V0);
    else
        require(false, "Invalid anode boundary condition: " + anodeBC);

    poisson_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.);
    poisson_solver.assemble_system_finalize();

    // solve Poisson equation
    return poisson_solver.solve();
}

template<int dim>
int PicSolver<dim>::update_positions() {
    for  (size_t i = 0; i < electrons.size(); ++i) {
        SuperParticle &particle = electrons.parts[i];

        //update position
        particle.pos += particle.vel * dt;

        //apply periodic boundaries
        particle.pos.x = periodic_image(particle.pos.x, box.xmax, box.xmin);
        particle.pos.y = periodic_image(particle.pos.y, box.ymax, box.ymin);

        // Update the cell ID; if any particles have left the domain their ID is set to -1
        // and they will be removed once we call clear_lost
        particle.cell = interpolator.update_point_cell(particle.pos, particle.cell);
    }

    int n_lost_particles = electrons.clear_lost();
    electrons.sort();
    return n_lost_particles;
}

template<int dim>
void PicSolver<dim>::update_velocities(){
    //update field
    for (auto particle : electrons.parts) {
        dealii::Point<dim> p(particle.pos.x, particle.pos.y, particle.pos.z);
        dealii::Tensor<1,dim> Efield = poisson_solver.probe_efield(p, particle.cell) ;

        //update velocities (corresponds to t + .5dt)
        particle.vel += Vec3(Efield) * (dt * electrons.q_over_m_factor) ;
    }
}

template<int dim>
int PicSolver<dim>::run_cycle(bool first_time) {
    // calculate electric field
    int n_cg_steps = compute_field(first_time);

    // updating the PIC particle velocities
    update_velocities();

    // collide PIC particles
    if (coll_coulomb_ee)
        this->coll_el_knm_2D();

    return n_cg_steps;
}

#define RAND ( (double)std::rand()/ RAND_MAX )
#define SQU(x)     ((x)*(x))
#define  PI     3.1415926535897932
#define  TWOPI  6.2831853071795864

template<int dim>
void PicSolver<dim>::coll_el_knm_2D() {
    ParticleSpecies &pa = electrons;

    static std::vector<size_t> inds2coll; //Not nice for parallelization
                                          // (neither over spieces- or cell)

    static const double Amplcoulomb =  1.;  // Amplification of coulomb collisions, only for testing purposes
    static const double LanLog = 13.;       // Coulomb Log
    //Constant factor in coulomb collisions
    //double Acoll = Amplcoulomb * ( !kind ? 1.0 : dt_ion*dt_ion*dt_ion ) * LanLog * SQU(SQU(Omega_pe)) * ncoll /
    //    ( TWOPI*PI * SQU(SQU(dz))*SQU(dz) * SQU(pa->mass) * SQU(Ndb) * N_sp );
    double Acoll = 1.0; //TODO !!!!!

    size_t Next = 0; // Index of the first particle in the current cell
    // (the way this is done today messes up paralellization)
    for (auto ord : pa.ordcount) {
        //Number of particles left to collide in the current cell
        size_t N2coll = ord;
        //Resize the indices
        inds2coll.resize(N2coll);

        //Constant factor in columb collisions for this cell
        double Acoll_cell = Acoll * ord / poisson_solver.get_cell_vol(pa.parts[Next].cell);

        if ( N2coll>1 ) {
            //Pick two random particle indices (j,k) from the particles that have not yet collide
            for (size_t loopind  = 0; loopind < N2coll; loopind++ ) {
                inds2coll[loopind]=loopind;
            }
            size_t j(0), k(0);
            while ( N2coll>0 ) {
                if ( N2coll>1 ) {
                    size_t jind = (size_t)(N2coll*RAND);
                    if (jind == N2coll) jind--;
                    j = inds2coll[jind];
                    for (size_t loopind = jind + 1; loopind < N2coll; loopind++) {
                        inds2coll[loopind-1] = inds2coll[loopind];
                    }
                    N2coll--;

                    jind = (size_t)(N2coll*RAND);
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
                Vec3 v_rel = pa.parts[j+Next].vel - pa.parts[k+Next].vel;
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
                double delta = (RAND-0.5)*deltaVar2;            // TODO !!!!!! Temporary hack until we have GausRandom
                double sinTheta = 2*delta/(1+SQU(delta));
                double oneMinusCosTheta = sinTheta*delta;
                //Scattering azimuth angle
                double Phi = RAND*TWOPI;
                double cosPhi = cos(Phi);
                double sinPhi = sin(Phi);

                //v_rel projected into XY plane (and inverse)
                Vec3 v_rel_DELTA;
                double v_rel_XY  = sqrt( SQU(v_rel.x) + SQU(v_rel.y) );
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
                pa.parts[j + Next].vel += v_rel_DELTA/2.;
                pa.parts[k + Next].vel -= v_rel_DELTA/2.;
            }//END do-loop over particles
        }//END if (N2coll > 1)
        Next += ord;
    }//END loop over cells
}
    
template<int dim>
void PicSolver<dim>::write(const string filename, const double time) const {
    string ftype = get_file_type(filename);
    require(ftype == "xyz" || ftype == "movie", "Invalid file type: " + ftype);

    if (electrons.size() == 0) return;  // Ovito can't handle 0 particles

    ofstream out;

    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    out << electrons.size() << endl;
    out << "time=" << time << ", PicSolver properties=id:I:1:pos:R:3:Velocity:R:3:cell:I:1" << endl;

    out.setf(ios::scientific);
    out.precision(6);

    for (int i = 0; i < electrons.size();  ++i)
        out << i << " " << electrons.parts[i] << endl;

    out.close();
}

/* ==================================================================
 *  ============================ PIC ===============================
 * ================================================================== */

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
int Pic<dim>::compute_field(bool first_time) {

    //TODO : save system lhs and copy it . don't rebuild it every time

    // Set-up system lhs
    laplace_solver.setup_system(first_time);
    if (first_time){
        laplace_solver.assemble_system_lhs();
        laplace_solver.save_system();
    } else
        laplace_solver.restore_system();

    // Set-up system RHS and BCs
    laplace_solver.assemble_system_pointcharge(electrons);
    if (anodeBC == "neumann")
        laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top, -E0);
    else if(anodeBC == "dirichlet")
        laplace_solver.assemble_system_dirichlet(fch::BoundaryId::vacuum_top, V0);
    else
        require(false, "Invalid anode boundary condition: " + anodeBC);

    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.);
    laplace_solver.assemble_system_finalize();

    // solve Poisson equation
    return laplace_solver.solve();
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
        particle.cell = interpolator.update_point_cell(particle.pos, particle.cell);
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
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(p, particle.cell) ;

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
        coll_el_knm_2D(electrons, dt, laplace_solver);

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
//template class Pic<2>;
template class Pic<3>;
template class PicSolver<3>;

} // namespace femocs
