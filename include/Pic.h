/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#ifndef PIC_H_
#define PIC_H_

#include "laplace.h"

#include <deal.II/base/point.h>

namespace femocs {


  template<int dim> class Pic {
  public:
    Pic();
    ~Pic();

    //Injects electrons
    // Indexing: (x1 y1 [z1] x2 y2 [z2] ...)
    int injectElectrons(const double* const r, const size_t n);
    
    //Computes the charge density for each FEM DOF
    int computeDensity();
    
    //Pushes the particles given the fields
    // - dt[s]
    void pushParticles(const double dt);
    
  private:

    //ELECTRONS
    //Particle positions [Å]
    std::vector<dealii::Point<dim>> r_el;
    //Particle velocities [Å/fs = 10 um/s]
    std::vector<dealii::Point<dim>> v_el;
    //Management
    std::vector<int> cid_el; //Index of the cell where the particle is inside

    //Constants
    const double q_over_m = 1.0; // [?] charge/mass for electrons
    const double q = 1.0; // [?] Charge of the particles (positive)
    
  };

}

#endif
