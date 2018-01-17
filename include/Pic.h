/*
 * Pic.h
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#ifndef PIC_H_
#define PIC_H_

#include "laplace.h"

namespace femocs {

  class Pic {
  public:
    Pic();
    ~Pic();

    //Injects electrons
    int injectElectrons(double* x, double* y, double* z, size_t n);
    
    //Computes the charge density for each FEM DOF
    int computeDensity();
    
    //Pushes the particles given the fields
    // dt[s]
    int pushParticles(const double dt);
    
  private:

    //ELECTRONS
    //Particle positions [Å]
    std::vector<double> x_el;
    std::vector<double> y_el;
    std::vector<double> z_el;
    //Particle velocities [Å/fs]
    std::vector<double> vx_el;
    std::vector<double> vy_el;
    std::vector<double> vz_el;
    //Management
    std::vector<int> cid_el; //Index of the cell where the particle is inside
  };

}

#endif
