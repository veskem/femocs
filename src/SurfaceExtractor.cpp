/*
 * SurfaceExtracter.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "SurfaceExtractor.h"

#include <omp.h>
#include <cmath>
#include <vector>

using namespace std;
namespace femocs {

// SurfaceExtractor constructor
SurfaceExtractor::SurfaceExtractor(const string& adapter, const double cutoff,
        const double latconst, const int nnn, const int nt) {
    this->adapter = adapter;
    this->cutoff2 = cutoff * cutoff;
    this->latconst = latconst;
    this->nnn = nnn;
    this->surfEps = cutoff / 1.0;
    this->nrOfOmpThreads = nt;
}

// Function to pick suitable extraction function
// Shared pointer prevents coping of data and makes memory management automatic
const shared_ptr<Surface> SurfaceExtractor::extract_surface(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {
    if (adapter == "coordination" && data->simu_type == "md")
        return coordination_extract(data, cell);
    else if (data->simu_type == "kmc")
        return kmc_extract(data, cell);
    else
        return shared_ptr<Surface>();
}

// Detect the atoms on the edges of horizontal plane or at the bottom of simulation cell
const bool SurfaceExtractor::is_boundary(const int i, const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {
    return (cell->xmax - abs(dat->x[i]) < surfEps) || (cell->ymax - abs(dat->y[i]) < surfEps)
            || (dat->z[i] - cell->zmin < surfEps);
}

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::kmc_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {

    int natoms = dat->x.size();
    int nsurface_atoms = 0;
    vector<bool> is_surface(natoms);
    int i;

    // Get number of surface atoms
    // cutoff == 0 turns off the surface extraction
    if (cutoff2 == 0) {
        for (i = 0; i < natoms; ++i)
            if (dat->type[i] != cell->type_vacancy) {
                is_surface[i] = 1;
                nsurface_atoms++;
            }
    } else {
        for (i = 0; i < natoms; ++i)
            if (dat->type[i] == cell->type_surf) {
                is_surface[i] = 1;
                nsurface_atoms++;
            }
    }
    // Create empty surface and add atoms
    // NOTICE, THAT SURFACE ATOM COORDINATION IS SET TO 0
    shared_ptr<Surface> s(new Surface(nsurface_atoms, this->latconst));
    for (i = 0; i < natoms; ++i)
        if (is_surface[i]) s->add(dat->x[i], dat->y[i], dat->z[i], 0, dat->type[i]);

    return s;
}

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::coordination_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {
    int natoms = dat->x.size();
    int i, j, coord;
    double r2, dx, dy, dz;
    vector<bool> is_surface(natoms);
    vector<int> coords(natoms);

    omp_set_num_threads(this->nrOfOmpThreads);

    int nsurface_atoms = 0;
    coord = 0;

    // cutoff == 0 turns off the surface extraction
    if (cutoff2 == 0) {
        for (i = 0; i < natoms; ++i)
            if (dat->type[i] != 3) {
                is_surface[i] = 1;
                coords[i] = 0;
                nsurface_atoms++;
            }

    } else {

#pragma omp parallel for shared(coords,is_surface) private(i,j,dx,dy,dz,r2) reduction(+:coord,nsurface_atoms)
        for (i = 0; i < natoms; ++i) {
            coord = 0;
            for (j = 0; j < natoms; ++j) {
                if (i == j) continue;
                dx = dat->x[i] - dat->x[j];
                dy = dat->y[i] - dat->y[j];
                dz = dat->z[i] - dat->z[j];
                r2 = dx * dx + dy * dy + dz * dz;
                if (r2 <= cutoff2) coord++;
                if (coord >= this->nnn) break; // Coordination can't be bigger than the biggest expected
            }
            coords[i] = coord;
            if ((coord < this->nnn) && (!is_boundary(i, dat, cell))) {
                nsurface_atoms++;
                is_surface[i] = 1;
            }
        }
    }

    // Create empty surface and add atoms by their coord
    shared_ptr<Surface> s(new Surface(nsurface_atoms, this->latconst));
    for (i = 0; i < natoms; ++i)
        if (is_surface[i]) s->add(dat->x[i], dat->y[i], dat->z[i], coords[i], dat->type[i]);

    return s;
}

// Extract surface by centrosymmetry analysis
const shared_ptr<Surface> SurfaceExtractor::centrosymmetry_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {
    return NULL;
}

} /* namespace femocs */
