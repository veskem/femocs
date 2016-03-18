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
    this->latconst = conf->latconst;
    this->nnn = conf->nnn;
    this->nrOfOmpThreads = conf->nt;
}


// Determine whether an atom is near the edge of simulation box
const bool SurfaceExtractor::on_edge(const double x, const double x_boundary) {
    return fabs(x - x_boundary) <= this->latconst / 2.3;
}

// Exctract the atoms near the simulation box sides
const shared_ptr<Surface> SurfaceExtractor::extract_edge(const shared_ptr<Surface> atoms,
        const Femocs::SimuCell* cell) {
    int N = atoms->get_n_atoms();
    int i, n_atoms;

    shared_ptr<Surface> edge(new Surface(N, this->latconst));

    for (i = 0; i < N; ++i) {
        if (on_edge(atoms->get_x(i), cell->xmin))
            edge->add_atom(cell->xmin, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i), cell->type_edge);
        if (on_edge(atoms->get_x(i), cell->xmax))
            edge->add_atom(cell->xmax, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i), cell->type_edge);
        if (on_edge(atoms->get_y(i), cell->ymin))
            edge->add_atom(atoms->get_x(i), cell->ymin, atoms->get_z(i), atoms->get_coordination(i), cell->type_edge);
        if (on_edge(atoms->get_y(i), cell->ymax))
            edge->add_atom(atoms->get_x(i), cell->ymax, atoms->get_z(i), atoms->get_coordination(i), cell->type_edge);
    }

    return edge;
}

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const shared_ptr<Surface> SurfaceExtractor::extract_reduced_bulk(const shared_ptr<Surface> surf,
        const Femocs::SimuCell* cell) {

    int i;
    auto edge = extract_edge(surf, cell);

    int nsurf = surf->get_n_atoms();
    int nedge = 4; //edge->getN();

    // Create empty surface and add atoms by their coord
    shared_ptr<Surface> bulk(new Surface(nsurf + nedge, this->latconst));
    for (i = 0; i < nsurf; ++i)
        bulk->add_atom(surf->get_x(i), surf->get_y(i), surf->get_z(i), surf->get_coordination(i), cell->type_surf);

//    for (i = 0; i < nedge; ++i)
//        bulk->add_atom(edge->getX(i), edge->getY(i), cell->zmin, 0, cell->type_surf);

    bulk->add_atom(cell->xmin, cell->ymin, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmin, cell->ymax, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmax, cell->ymin, cell->zmin, 0, cell->type_bulk);
    bulk->add_atom(cell->xmax, cell->ymax, cell->zmin, 0, cell->type_bulk);

    bulk->calc_statistics();
    return bulk;
}

// Function to extract bulk material from input atomistic data
const shared_ptr<Surface> SurfaceExtractor::extract_truncated_bulk(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {

    int i;
    int N = data->x.size();
    vector<bool> is_bulk(N);
    vector<bool> is_surf(N);

    // Get number and locations of bulk atoms
    for (i = 0; i < N; ++i)
        is_bulk[i] = (data->type[i] == cell->type_bulk || data->type[i] == cell->type_surf)
            && (data->z[i] > cell->zmin);
        //is_bulk[i] = (data->z[i] > cell->zmin) && (data->type[i] == cell->type_bulk);

    for (i = 0; i < N; ++i)
        is_surf[i] = (data->coordination[i] > 0) && (data->coordination[i] < this->nnn);

    int M = accumulate(is_bulk.begin(), is_bulk.end(), 0);
    shared_ptr<Surface> bulk(new Surface(M, this->latconst));

    // Add bulk atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && !is_surf[i])
            bulk->add_atom(data->x[i], data->y[i], data->z[i], data->coordination[i], cell->type_bulk);

    // Add surface atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && is_surf[i])
            bulk->add_atom(data->x[i], data->y[i], data->z[i], data->coordination[i], cell->type_surf);

    bulk->calc_statistics();
    return bulk;
}

// Function to extract bulk material from input atomistic data
const shared_ptr<Surface> SurfaceExtractor::extract_bulk(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {
    int i;
    int N = data->x.size();
    vector<bool> is_bulk(N);

    for (i = 0; i < N; ++i)
        is_bulk[i] = data->type[i] != cell->type_vacancy;

    int M = accumulate(is_bulk.begin(), is_bulk.end(), 0);
    shared_ptr<Surface> bulk(new Surface(M, this->latconst));
    for (i = 0; i < N; ++i)
        if (is_bulk[i])
            bulk->add_atom(data->x[i], data->y[i], data->z[i], data->coordination[i], cell->type_bulk);

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

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::kmc_extract(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {

    int i;
    int N = data->x.size();
    vector<bool> is_surface(N);

    // Get number and locations of surface atoms
    for (i = 0; i < N; ++i)
        is_surface[i] = (data->type[i] == cell->type_surf);

    int M = accumulate(is_surface.begin(), is_surface.end(), 0);
    // Create empty Surface
    shared_ptr<Surface> surf(new Surface(M, this->latconst));

    // Add surface atoms to Surface
    for (i = 0; i < N; ++i)
        if (is_surface[i])
            surf->add_atom(data->x[i], data->y[i], data->z[i], data->coordination[i], cell->type_surf);

    surf->calc_statistics();
    return surf;
}

// Extract surface by coordination analysis
const shared_ptr<Surface> SurfaceExtractor::coordination_extract(const AtomReader::Data* data,
        const Femocs::SimuCell* cell) {

    int N = data->x.size();
    int i, j, coord;
    double r2, dx, dy, dz;
    double zmin = cell->zmin + this->latconst;

    vector<bool> is_surface(N);
    vector<int> coords(N);

    for (i = 0; i < N; ++i)
        is_surface[i] = (data->coordination[i] > 0) && (data->coordination[i] < this->nnn)
                && (data->z[i] > zmin);

    // Create empty surface and add atoms by their coord
    int M = accumulate(is_surface.begin(), is_surface.end(), 0);
    shared_ptr<Surface> surf(new Surface(M, this->latconst));

    for (i = 0; i < N; ++i)
        if (is_surface[i])
            surf->add_atom(data->x[i], data->y[i], data->z[i], data->coordination[i],
                    cell->type_surf);

    surf->calc_statistics();
    return surf;
}

// Function to extract bulk material from input atomistic data
const void SurfaceExtractor::rectangularize(shared_ptr<Surface> atoms,
        const Femocs::SimuCell* cell) {

    double zmin_up = atoms->sizes.zmin + atoms->sizes.latconst/2.3;
    double zmin_down = atoms->sizes.zmin;

    for (int i = 0; i < atoms->get_n_atoms(); ++i) {
        // Flatten the atoms on bottom layer
        if( (atoms->get_z(i) >= zmin_down) && (atoms->get_z(i) <= zmin_up) ) atoms->set_z(i, zmin_down);

        // Flatten the atoms on the sides of simulation box
        if (on_edge(atoms->get_x(i), cell->xmin)) atoms->set_x(i, cell->xmin);
        if (on_edge(atoms->get_x(i), cell->xmax)) atoms->set_x(i, cell->xmax);
        if (on_edge(atoms->get_y(i), cell->ymin)) atoms->set_y(i, cell->ymin);
        if (on_edge(atoms->get_y(i), cell->ymax)) atoms->set_y(i, cell->ymax);
    }

    // Add atoms to the bottom corner of the simulation cell
//    atoms->add_atom(atoms->sizes.xmin, atoms->sizes.ymin, zmin_down, 0, cell->type_bulk);
//    atoms->add_atom(atoms->sizes.xmin, atoms->sizes.ymax, zmin_down, 0, cell->type_bulk);
//    atoms->add_atom(atoms->sizes.xmax, atoms->sizes.ymin, zmin_down, 0, cell->type_bulk);
//    atoms->add_atom(atoms->sizes.xmax, atoms->sizes.ymax, zmin_down, 0, cell->type_bulk);
}

} /* namespace femocs */
