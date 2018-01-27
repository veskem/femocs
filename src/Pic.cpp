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
template<int dim>
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver, FieldReader &fr,
        fch::CurrentsAndHeating<3> &ch_solver, HeatReader &hr, EmissionReader &er) : ///< Object to read the temperature data) :
        laplace_solver(laplace_solver), fr(fr), ch_solver(ch_solver), hr(hr), er(er){}

template<int dim>
Pic<dim>::~Pic() {}

template<int dim>
int Pic<dim>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        r_el.push_back(dealii::Point<3>(r[i*3+0],r[i*3+1],r[i*3+2]));
        v_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
        F_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
        cid_el.push_back(fr.update_point_cell(r_el[i], 10));
    }
}

template<int dim>
int Pic<dim>::injectElectrons() {
    vector<dealii::Point<dim>> positions, fields;
    vector<int> cells;
    er.inject_electrons(dt, Wsp, positions, fields, cells);

    for (int i = 0; i < fields.size(); ++i){
        r_el.push_back(positions[i]);
        F_el.push_back(-fields[i]);
        v_el.push_back(F_el.back() * q_over_m_factor * dt * .5);
        cid_el.push_back(cells[i]);
    }
    return 0;

}

//Call the laplace solver with the list of positions and charge(s)
template<int dim>
void Pic<dim>::computeField() {

    double t0;
    start_msg(t0, "=== Solving the Poisson equation...");
    laplace_solver.setup_system();
    laplace_solver.assemble_system_lhs();
    bool boundary_set = false;

    if (anodeBC == "neumann"){
        cout << "assembling neumann boundary for E0 = " << E0 << endl;
        laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top, -E0);
        boundary_set = true;
    }

    laplace_solver.assemble_system_pointcharge(r_el, -q_over_eps0*Wsp, cid_el);

    if(anodeBC == "dirichlet"){
        cout << "assembling dirichlet boundary" << endl;
        laplace_solver.assemble_system_dirichlet(fch::BoundaryId::vacuum_top, V0);
        boundary_set = true;
    }

    require(boundary_set, "ERROR: anodeBC parameter wrong!! anodeBC = " + anodeBC);

    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.);

    laplace_solver.assemble_system_finalize();
    laplace_solver.solve();
    end_msg(t0);
}


template<int dim>
void Pic<dim>::updatePositions(){
    for (size_t i = 0; i < r_el.size(); i++) {

        //update position
        r_el[i] = r_el[i] + v_el[i]*dt;

        //Update the cid_el && check if any particles have left the domain && remove them
        cid_el[i] = fr.update_point_cell(r_el[i], cid_el[i]);
        if (cid_el[i] == -1) {
            lost_el.push_back(i);
        }
    }
}

template<int dim>
void Pic<dim>::updateFieldAndVelocities(){

    //update field
    for (size_t i = 0; i < v_el.size(); i++) {
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(r_el[i], cid_el[i]) ;
        F_el[i] = -dealii::Point<dim>(Efield) ;

        //update velocities (corresponds to t + .5dt)
        v_el[i] = v_el[i] + q_over_m_factor*F_el[i]*(dt);

    }
}

template<int dim>
void Pic<dim>::runCycle() {
    updatePositions();
    clearLostParticles();
    computeField();
    updateFieldAndVelocities();
}

template<int dim>
void Pic<dim>::clearLostParticles(){
    size_t npart = r_el.size();
    size_t nlost = 0;

    //Delete the lost particles from the arrays
    for (size_t i = 0; i < npart; i++) {
        bool islost=false;
        //Is this particle lost?
        for (auto lost : lost_el) {
            if (lost == i) {
                islost=true;
                nlost++;
                break;
            }
        }
        if (nlost==0 or islost) continue; // Don't shuffle this particle left

        r_el[i-nlost] = r_el[i];
        v_el[i-nlost] = v_el[i];
        F_el[i-nlost] = F_el[i];
        cid_el[i-nlost] = cid_el[i];
    }

    //Shrink the arrays
    if (nlost > 0){
        r_el.resize(npart-nlost);
        v_el.resize(npart-nlost);
        cid_el.resize(npart-nlost);
        F_el.resize(npart-nlost);
        cout << "Particles where lost! nlost=" << nlost << endl;
    }

    lost_el.clear();
}

template<int dim>
void Pic<dim>::writeParticles(const string filename) {

    // dummy electron always added in the end (to avoid empty electrons crashing ovito)
    ofstream out;
    out.setf(std::ios::scientific);
    out.precision(6);

    string ftype = get_file_type(filename);
    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    cout << "writing particles to " + filename << " n_size = " << r_el.size() << endl;

    out << r_el.size() + 1 << endl;
    out << "Interpolator properties=id:I:1:pos:R:3:vel:R:3:Force:R:3:cell:I:1" << endl;

    for (int i = 0; i < r_el.size(); ++i)
        out << i << " " << r_el[i][0] << " " << r_el[i][1] << " " << r_el[i][2] << " " <<
        v_el[i][0] << " " << v_el[i][1] << " " << v_el[i][2] << " " <<
        F_el[i][0] << " " << F_el[i][1] << " " << F_el[i][2] << " " << cid_el[i] << endl;

    out << "-1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0" << endl;

    out.close();
}


//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;
}
