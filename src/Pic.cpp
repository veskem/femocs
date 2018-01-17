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

  template<int dim> int Pic<dim>::injectElectrons(const double* const r, const size_t n) {
    
  }
  
  template<int dim> int Pic<dim>::computeDensity() {

  }

  template<int dim> int Pic<dim>::pushParticles(const double dt) {
    
  }

  //Tell the compiler which types to actually compile, so that they are available for the linker
  template class Pic<2>;
  template class Pic<3>;
}

