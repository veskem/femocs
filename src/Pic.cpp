/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#include "Pic.h"

#include <deal.II/base/tensor.h>

namespace femocs {
template<int dim>
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver) : laplace_solver(laplace_solver){
}

template<int dim>
Pic<dim>::~Pic() {

}

template<>
int Pic<2>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        r_el.push_back(dealii::Point<2>(r[i*2+0],r[i*2+1]));
        v_el.push_back(dealii::Point<2>(0.0,0.0));
        cid_el.push_back(-1);//Unknown cell
    }
}
template<>
int Pic<3>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        r_el.push_back(dealii::Point<3>(r[i*3+0],r[i*3+1],r[i*3+2]));
        v_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
        cid_el.push_back(-1);//Unknown cell
    }
}

template<int dim>
void Pic<dim>::computeField() {
    //Call the laplace solver with the list of positions and charge(s)
    laplace_solver.assemble_system_lhs();
    laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top);

    laplace_solver.assemble_system_pointcharge(r_el, q, cid_el);
    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.0);
}

template<int dim>
void Pic<dim>::pushParticles(const double dt, FieldReader &fr) {


    for (size_t i = 0; i < r_el.size(); i++) {
        //Leapfrog method:
        // positions defined ON the time steps, velocities defined at half time steps
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(r_el[i], cid_el[i]) ; // Get the field!

        v_el[i] = v_el[i] + q_over_m*Efield*dt;
        r_el[i] = r_el[i] + v_el[i]*dt;

        //Update the cid_el && check if any particles have left the domain
        cid_el[i] = fr.update_point_cell(r_el[i], cid_el[i]);
    }

}

//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;
}

