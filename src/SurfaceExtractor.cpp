/*
 * SurfaceExtracter.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "SurfaceExtractor.h"

#include <omp.h>
#include <cmath>
#include <numeric>
#include <vector>

using namespace std;
namespace femocs {

// SurfaceExtractor constructor
SurfaceExtractor::SurfaceExtractor(const Femocs::Config* conf) {
    this->adapter = conf->extracter;
    this->cutoff2 = conf->coord_cutoff * conf->coord_cutoff;
    this->latconst = conf->latconst;
    this->nnn = conf->nnn;
    this->surfEps = conf->coord_cutoff / 1.0;
    this->nrOfOmpThreads = conf->nt;
}

const shared_ptr<Surface> SurfaceExtractor::extract_edge(const shared_ptr<Surface> surf, const Femocs::SimuCell* cell) {
    int N = surf->getN();
    int i, natoms;
    vector<bool> is_edge(N);
    is_edge.reserve(N);

    for (i = 0; i < N; ++i)
        is_edge[i] = on_edge(surf->getX(i), surf->getY(i), cell);

    natoms = accumulate(is_edge.begin(), is_edge.end(), 0);
    shared_ptr<Surface> edge(new Surface(natoms, this->latconst));

    for (i = 0; i < N; ++i)
        if (is_edge[i])
            edge->add(surf->getX(i), surf->getY(i), surf->getZ(i), 0, cell->type_edge);

    return edge;
}

// !!! REPLACE WITH COORDINATION DATA !!!
const bool SurfaceExtractor::on_edge(const double x, const double y, const Femocs::SimuCell* cell) {
    const double eps = 0.1;
    bool dif1 = fabs(x - cell->xmin) < eps || fabs(x - cell->xmax) < eps;
    bool dif2 = fabs(y - cell->ymin) < eps || fabs(y - cell->ymax) < eps;
    return dif1 | dif2;
}

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const shared_ptr<Surface> SurfaceExtractor::extract_bulk_reduced(const shared_ptr<Surface> surf,
        const Femocs::SimuCell* cell) {

    int i;
    auto edge = extract_edge(surf, cell);

    int nsurf = surf->getN();
    int nedge = edge->getN();

    // Create empty surface and add atoms by their coord
    shared_ptr<Surface> bulk(new Surface(nsurf+nedge, this->latconst));
    for (i = 0; i < nsurf; ++i)
        bulk->add(surf->getX(i), surf->getY(i), surf->getZ(i), 0, cell->type_surf);

    for (i = 0; i < nedge; ++i)
        bulk->add(edge->getX(i), edge->getY(i), cell->xmin, 0, cell->type_surf);

    return bulk;
}

// Function to extract bulk material from input atomistic data
const shared_ptr<Surface> SurfaceExtractor::extract_bulk(const AtomReader::Data* data, const Femocs::SimuCell* cell) {
    int i;
    int N = data->x.size();
    vector<bool> is_quality(N);

    for (i = 0; i < N; ++i)
        is_quality[i] = data->type[i] != cell->type_vacancy;

    int M = accumulate(is_quality.begin(), is_quality.end(), 0);
    shared_ptr<Surface> bulk( new Surface(M, this->latconst) );
    for (i = 0; i < N; ++i)
        if ( is_quality[i] )
            bulk->add(data->x[i], data->y[i], data->z[i], 0, cell->type_bulk);

    return bulk;
}

// Function to pick suitable extraction function
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

    int N = dat->x.size();
    vector<bool> is_surface(N);
    int i;

    // Get number and locations of surface atoms
    for (i = 0; i < N; ++i)
        is_surface[i] = (dat->type[i] == cell->type_surf);

    int M = accumulate(is_surface.begin(), is_surface.end(), 0);
    // Create empty surface and add atoms
    shared_ptr<Surface> surf(new Surface(M, this->latconst));
    for (i = 0; i < N; ++i)
        if (is_surface[i])
            surf->add(dat->x[i], dat->y[i], dat->z[i], 0, dat->type[i]);

    return surf;
}

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::coordination_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {

    int N = dat->x.size();
    int i, j, coord;
    double r2, dx, dy, dz;
    vector<bool> is_surface(N);
    vector<int> coords(N);

    omp_set_num_threads(this->nrOfOmpThreads);

    coord = 0;

#pragma omp parallel for shared(coords,is_surface) private(i,j,dx,dy,dz,r2) reduction(+:coord,N)
    for (i = 0; i < N; ++i) {
        coord = 0;
        for (j = 0; j < N; ++j) {
            if (i == j) continue;
            dx = dat->x[i] - dat->x[j];
            dy = dat->y[i] - dat->y[j];
            dz = dat->z[i] - dat->z[j];
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 <= cutoff2) coord++;
            if (coord >= this->nnn) break; // Coordination can't be bigger than the biggest expected
        }
        coords[i] = coord;
        is_surface[i] = ( (coord < this->nnn) && (!is_boundary(i, dat, cell)) );
    }

    // Create empty surface and add atoms by their coord
    int M = accumulate(is_surface.begin(), is_surface.end(), 0);
    shared_ptr<Surface> surf(new Surface(M, this->latconst));
    for (i = 0; i < N; ++i)
        if ( is_surface[i] )
            surf->add(dat->x[i], dat->y[i], dat->z[i], coords[i], dat->type[i]);

    return surf;
}

// Extract surface by centrosymmetry analysis
const shared_ptr<Surface> SurfaceExtractor::centrosymmetry_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {
    return NULL;
}

} /* namespace femocs */
