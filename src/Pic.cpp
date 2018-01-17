/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#include "Pic.h"

namespace femocs {
  template<int dim> Pic<dim>::Pic() {

  }

  template<int dim> Pic<dim>::~Pic() {

  }

  template<> int Pic<2>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
      r_el.push_back(dealii::Point<2>(r[i*2+0],r[i*2+1]));
      v_el.push_back(dealii::Point<2>(0.0,0.0));
      cid_el.push_back(-1);//Unknown cell
    }
  }
  template<> int Pic<3>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
      r_el.push_back(dealii::Point<3>(r[i*3+0],r[i*3+1],r[i*3+2]));
      v_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
      cid_el.push_back(-1);//Unknown cell
    }
  }
  
  template<int dim> int Pic<dim>::computeDensity() {
    //Call the laplace solver with the list of positions and charge(s)
  }

  template<int dim> void Pic<dim>::pushParticles(const double dt) {

    for (size_t i = 0; i < r_el.size(); i++) {
      //Leapfrog method:
      // positions defined ON the time steps, velocities defined at half time steps
      dealii::Point<dim> Efield; // Get the field!
      
      v_el[i] = v_el[i] + q_over_m*Efield*dt;
      r_el[i] = r_el[i] + v_el[i]*dt;

      //Update the cid_el && check if any particles have left the domain
      
    }

  }

  //Tell the compiler which types to actually compile, so that they are available for the linker
  template class Pic<2>;
  template class Pic<3>;
}

