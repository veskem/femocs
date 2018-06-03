/*
 * PicCollisions.cpp
 *
 *  Created on: 31.01.2018
 *      Author: Kyrre
 */

#include "PicCollisions.h"

#define RAND ( (double)std::rand()/ RAND_MAX )

#define SQU(x)     ((x)*(x))
#define  PI     3.1415926535897932
#define  TWOPI  6.2831853071795864

namespace femocs {
    void coll_el_knm_2D( ParticleSpecies &pa, const double dt, PoissonSolver<3> &poisson_solver ) {

        static std::vector<size_t> inds2coll; //Not nice for parallelization
                                              // (neither over spieces- or cell)

//        static const double Amplcoulomb =  1.;  // Amplification of coulomb collisions, only for testing purposes
//        static const double LanLog = 13.;       // Coulomb Log
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
}
