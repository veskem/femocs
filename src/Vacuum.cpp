/*
 * Vacuum.cpp
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#include "Vacuum.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;
namespace femocs {

const void Vacuum::initialise_vars(const int natoms) {
    this->x.reserve(natoms);
    this->y.reserve(natoms);
    this->z.reserve(natoms);
    this->type.reserve(natoms);
    this->N = 0;
}

Vacuum::Vacuum() {
}

const void Vacuum::generate_simple(const Femocs::SimuCell* cell) {
    int M = 4; // total number of nodes
    initialise_vars(M);

    // Add points to the xy-plane corners on current layer
    addpoint(cell->xmax, cell->ymax, cell->zmaxbox, cell->type_vacuum);
    addpoint(cell->xmax, cell->ymin, cell->zmaxbox, cell->type_vacuum);
    addpoint(cell->xmin, cell->ymax, cell->zmaxbox, cell->type_vacuum);
    addpoint(cell->xmin, cell->ymin, cell->zmaxbox, cell->type_vacuum);
}

const void Vacuum::generate(const Femocs::SimuCell* cell, shared_ptr<Surface> surf) {
    const double layer_exp = 1.7; // Magic number used to get nicer mesh distribution in z-direction
    const double zmax = cell->zmaxbox;
    const double dz = surf->getLatconst();
    int nlayers = ceil(log(zmax / dz) / log(layer_exp));  // number of node layers
    int N_surf = surf->getN();                          // number of atoms in surface

    // sum of Nsurf/2^1 + Nsurf/2^2 +...+ Nsurf/2^nlayers = Nsurf*( 1 - 1/2^nlayers )/(1 - 1/2)
    int M = ceil(4 * nlayers + N_surf * (1.0 - 1.0 / (1 << nlayers))); // total number of nodes
    initialise_vars(M);

    double znew = zmax;
    int i, j, k;

    vector<int> iperm;
    //iperm.reserve(N_surf);
    this->sort_atoms(&iperm);

    // Add every 2nd atom to the next layer as compared to previous one
    for (i = 1; i <= nlayers; ++i) {
        for (j = 0; j < N_surf / (1 << i); ++j) {

            // !!! SOME PROBLEMS WITH SORTING !!!
            //k = iperm[j * (1 << i)];
            k = j * (1 << i);

            if (k >= N_surf) {
                cout << "k out of Surface limits!" << endl;
                exit(EXIT_FAILURE);
            }

            znew = dz * pow(layer_exp, i) + surf->getZ(k);
            if (znew > zmax) znew = zmax;

            addpoint(surf->getX(k), surf->getY(k), znew, cell->type_vacuum);
        }

        // Add points to the xy-plane corners on current layer
        addpoint(cell->xmax, cell->ymax, znew, cell->type_vacuum);
        addpoint(cell->xmax, cell->ymin, znew, cell->type_vacuum);
        addpoint(cell->xmin, cell->ymax, znew, cell->type_vacuum);
        addpoint(cell->xmin, cell->ymin, znew, cell->type_vacuum);
    }
}

const void Vacuum::addpoint(const double x, const double y, const double z, const int type) {
//    if (this->N >= this->x.size()) {
//        cout << "Size of Vacuum exceeded!" << endl;
//    }
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->type.push_back(type);
    this->N++;
}

} /* namespace femocs */
