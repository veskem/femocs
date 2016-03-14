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

const shared_ptr<Surface> SurfaceExtractor::extract_edge(const shared_ptr<Surface> atoms,
        const Femocs::SimuCell* cell) {
    int N = atoms->get_N();
    int i, n_atoms;
    vector<bool> is_xmin(N);
    vector<bool> is_xmax(N);
    vector<bool> is_ymin(N);
    vector<bool> is_ymax(N);

    for (i = 0; i < N; ++i) {
        is_xmin[i] = on_edge_vol2(atoms->get_x(i), cell->xmin);
        is_xmax[i] = on_edge_vol2(atoms->get_x(i), cell->xmax);
        is_ymin[i] = on_edge_vol2(atoms->get_y(i), cell->ymin);
        is_ymax[i] = on_edge_vol2(atoms->get_y(i), cell->ymax);
    }

    n_atoms = accumulate(is_xmin.begin(), is_xmin.end(), 0)
            + accumulate(is_xmax.begin(), is_xmax.end(), 0)
            + accumulate(is_ymin.begin(), is_ymin.end(), 0)
            + accumulate(is_xmax.begin(), is_xmax.end(), 0);

    shared_ptr<Surface> edge(new Surface(n_atoms, this->latconst));

    for (i = 0; i < N; ++i) {
//        if (is_xmin[i]) edge->add(cell->xmin, atoms->getY(i), atoms->sizes.zmin, 0, cell->type_edge);
//        if (is_xmax[i]) edge->add(cell->xmax, atoms->getY(i), atoms->sizes.zmin, 0, cell->type_edge);
//        if (is_ymin[i]) edge->add(atoms->getX(i), cell->ymin, atoms->sizes.zmin, 0, cell->type_edge);
//        if (is_ymax[i]) edge->add(atoms->getX(i), cell->ymax, atoms->sizes.zmin, 0, cell->type_edge);

        if (is_xmin[i]) edge->add_atom(cell->xmin, atoms->get_y(i), atoms->get_z(i), 0, cell->type_edge);
        if (is_xmax[i]) edge->add_atom(cell->xmax, atoms->get_y(i), atoms->get_z(i), 0, cell->type_edge);
        if (is_ymin[i]) edge->add_atom(atoms->get_x(i), cell->ymin, atoms->get_z(i), 0, cell->type_edge);
        if (is_ymax[i]) edge->add_atom(atoms->get_x(i), cell->ymax, atoms->get_z(i), 0, cell->type_edge);
    }

    return edge;
}

// TODO !!! REPLACE WITH COORDINATION DATA !!!
// MAKE IT INLINE
const bool SurfaceExtractor::on_edge_vol2(const double x, const double x_boundary) {
    const double eps = 1.0;
    return fabs(x - x_boundary) <= eps;
}

// TODO !!! REPLACE WITH COORDINATION DATA !!!
const bool SurfaceExtractor::on_edge(const double x, const double y, const Femocs::SimuCell* cell) {
    const double eps = 0.1;
    bool dif1 = fabs(x - cell->xmin) < eps || fabs(x - cell->xmax) < eps;
    bool dif2 = fabs(y - cell->ymin) < eps || fabs(y - cell->ymax) < eps;
    return dif1 | dif2;
}

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const shared_ptr<Surface> SurfaceExtractor::extract_reduced_bulk(const shared_ptr<Surface> surf,
        const Femocs::SimuCell* cell) {

    int i;
    auto edge = extract_edge(surf, cell);

    int nsurf = surf->get_N();
    int nedge = 4; //edge->getN();

    // Create empty surface and add atoms by their coord
    shared_ptr<Surface> bulk(new Surface(nsurf + nedge, this->latconst));
    for (i = 0; i < nsurf; ++i)
        bulk->add_atom(surf->get_x(i), surf->get_y(i), surf->get_z(i), 0, cell->type_surf);

//    for (i = 0; i < nedge; ++i)
//        bulk->add(edge->getX(i), edge->getY(i), cell->zmin, 0, cell->type_surf);

    bulk->add_atom(cell->xmin, cell->ymin, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmin, cell->ymax, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmax, cell->ymin, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmax, cell->ymax, cell->zmin, 0, cell->type_bulk);

    bulk->calc_statistics();
    return bulk;
}

// Function to extract bulk material from input atomistic data
const shared_ptr<Surface> SurfaceExtractor::extract_truncated_bulk(const AtomReader::Data* data, const Femocs::SimuCell* cell) {

    int i;
    int N = data->x.size();
    vector<bool> is_bulk(N);

    for (i = 0; i < N; ++i)
        is_bulk[i] = (data->type[i] != cell->type_vacancy) && (data->z[i] > cell->zmin);

    int M = accumulate(is_bulk.begin(), is_bulk.end(), 0);

    shared_ptr<Surface> bulk(new Surface(M, this->latconst));

    // Add initial bulk atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i]) bulk->add_atom(data->x[i], data->y[i], data->z[i], 0, cell->type_bulk);

    bulk->calc_statistics();
    return bulk;
}

// Function to extract bulk material from input atomistic data
const void SurfaceExtractor::rectangularize_bulk(shared_ptr<Surface> bulk, const Femocs::SimuCell* cell) {

    double zmin_up   = bulk->sizes.zmin + bulk->sizes.latconst;
    double zmin_down = bulk->sizes.zmin;

    for (int i = 0; i < bulk->get_N(); ++i) {
        // Flatten the atoms on bottom layer
//        if( (bulk->get_z(i) >= zmin_down) && (bulk->get_z(i) <= zmin_up) ) bulk->set_z(i, bulk->sizes.zmin + 0.5*bulk->sizes.latconst);

        // Flatten the atoms on the sides of simulation box
        if( on_edge_vol2(bulk->get_x(i), cell->xmin) ) bulk->set_x(i, cell->xmin);
        if( on_edge_vol2(bulk->get_x(i), cell->xmax) ) bulk->set_x(i, cell->xmax);
        if( on_edge_vol2(bulk->get_y(i), cell->ymin) ) bulk->set_y(i, cell->ymin);
        if( on_edge_vol2(bulk->get_y(i), cell->ymax) ) bulk->set_y(i, cell->ymax);
    }

    // Add atoms to the bottom corner of the simulation cell
//    bulk->add(bulk->sizes.xmin, bulk->sizes.ymin, zmin_down, 0, cell->type_bulk);
//    bulk->add(bulk->sizes.xmin, bulk->sizes.ymax, zmin_down, 0, cell->type_bulk);
//    bulk->add(bulk->sizes.xmax, bulk->sizes.ymin, zmin_down, 0, cell->type_bulk);
//    bulk->add(bulk->sizes.xmax, bulk->sizes.ymax, zmin_down, 0, cell->type_bulk);

    //bulk->calc_statistics();
}

// Function to extract bulk material from input atomistic data
const shared_ptr<Surface> SurfaceExtractor::extract_bulk(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {
    int i;
    int N = data->x.size();
    vector<bool> is_quality(N);

    for (i = 0; i < N; ++i)
        is_quality[i] = data->type[i] != cell->type_vacancy;

    int M = accumulate(is_quality.begin(), is_quality.end(), 0);
    shared_ptr<Surface> bulk(new Surface(M, this->latconst));
    for (i = 0; i < N; ++i)
        if (is_quality[i]) bulk->add_atom(data->x[i], data->y[i], data->z[i], 0, cell->type_bulk);

    bulk->calc_statistics();
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
        if (is_surface[i]) surf->add_atom(dat->x[i], dat->y[i], dat->z[i], 0, dat->type[i]);

    surf->calc_statistics();
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
        is_surface[i] = ((coord < this->nnn) && (!is_boundary(i, dat, cell)));
    }

    // Create empty surface and add atoms by their coord
    int M = accumulate(is_surface.begin(), is_surface.end(), 0);
    shared_ptr<Surface> surf(new Surface(M, this->latconst));
    for (i = 0; i < N; ++i)
        if (is_surface[i]) surf->add_atom(dat->x[i], dat->y[i], dat->z[i], coords[i], dat->type[i]);

    surf->calc_statistics();
    return surf;
}

// Extract surface by centrosymmetry analysis
const shared_ptr<Surface> SurfaceExtractor::centrosymmetry_extract(const AtomReader::Data* dat,
        const Femocs::SimuCell* cell) {
    return NULL;
}

} /* namespace femocs */
