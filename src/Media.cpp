/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Media.h"

using namespace std;

namespace femocs {

// =================================================
//             Implementation of Vacuum
// =================================================

// Constructor of Vacuum class
Vacuum::Vacuum() {
    init_statistics();
}

// Generates Vacuum by adding four points to the top of simulation cell
const void Vacuum::generate_simple(const Femocs::SimuCell* cell) {
    int M = 4; // total number of nodes
    reserve(M);

    // Add points to the xy-plane corners on current layer
    add_atom(cell->xmax, cell->ymax, cell->zmaxbox, 0);
    add_atom(cell->xmax, cell->ymin, cell->zmaxbox, 0);
    add_atom(cell->xmin, cell->ymax, cell->zmaxbox, 0);
    add_atom(cell->xmin, cell->ymin, cell->zmaxbox, 0);

    init_statistics();
    calc_statistics();
}

// =================================================
//              Implementation of Bulk
// =================================================

// Constructor for Bulk class
Bulk::Bulk(const double latconst, const int nnn) {
    init_statistics();
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const void Bulk::extract_reduced_bulk(Surface* surf, const Femocs::SimuCell* cell) {
    int i;
    shared_ptr<Edge> edge(new Edge(crys_struct.latconst, crys_struct.nnn));
    edge->extract_edge(surf, cell);

    int n_surf = surf->get_n_atoms();
    int n_edge = 4; //edge->getN();

    // Reserve memory for bulk atoms and their parameters
    reserve(n_surf + n_edge);

    for (i = 0; i < n_surf; ++i)
        add_atom(surf->get_x(i), surf->get_y(i), surf->get_z(i), surf->get_coordination(i));

//    for (i = 0; i < n_edge; ++i)
//        add_atom(edge->get_x(i), edge->get_y(i), cell->zmin, edge->get_coordination(i));

    // Add extra atoms to the bottom corner edges of simulation cell
    add_atom(cell->xmin, cell->ymin, cell->zmin, 0);
    add_atom(cell->xmin, cell->ymax, cell->zmin, 0);
    add_atom(cell->xmax, cell->ymin, cell->zmin, 0);
    add_atom(cell->xmax, cell->ymax, cell->zmin, 0);

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract_truncated_bulk(AtomReader* reader, const Femocs::SimuCell* cell) {
    int i;
    int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);
    vector<bool> is_surf(N);

    // Get number and locations of bulk atoms
    for (i = 0; i < N; ++i)
        is_bulk[i] =
                (reader->get_z(i) > cell->zmin)
                        && (reader->get_type(i) == cell->type_bulk
                                || reader->get_type(i) == cell->type_surf);

    for (i = 0; i < N; ++i)
        is_surf[i] = (reader->get_coord(i) > 0) && (reader->get_coord(i) < crys_struct.nnn);

    reserve(accumulate(is_bulk.begin(), is_bulk.end(), 0));

    // Add bulk atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && !is_surf[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i), reader->get_coord(i));

    // Add surface atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && is_surf[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i), reader->get_coord(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract_bulk(AtomReader* reader, const Femocs::SimuCell* cell) {
    int i;
    int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);

    for (i = 0; i < N; ++i)
        is_bulk[i] = reader->get_type(i) != cell->type_vacancy;

    reserve(accumulate(is_bulk.begin(), is_bulk.end(), 0));

    for (i = 0; i < N; ++i)
        if (is_bulk[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i), reader->get_coord(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::rectangularize(const Femocs::SimuCell* cell) {
    double zmin_up = sizes.zmin + crys_struct.latconst / 2.1;
    double zmin_down = sizes.zmin;

    for (int i = 0; i < get_n_atoms(); ++i) {
        // Flatten the atoms on bottom layer
        if ((get_z(i) >= zmin_down) && (get_z(i) <= zmin_up)) set_z(i, zmin_down);

        // Flatten the atoms on the sides of simulation box
        if (on_edge(get_x(i), cell->xmin)) set_x(i, cell->xmin);
        if (on_edge(get_x(i), cell->xmax)) set_x(i, cell->xmax);
        if (on_edge(get_y(i), cell->ymin)) set_y(i, cell->ymin);
        if (on_edge(get_y(i), cell->ymax)) set_y(i, cell->ymax);
    }

    // Add atoms to the bottom corner of the simulation cell
//    add_atom(sizes.xmin, sizes.ymin, zmin_down, 0);
//    add_atom(sizes.xmin, sizes.ymax, zmin_down, 0);
//    add_atom(sizes.xmax, sizes.ymin, zmin_down, 0);
//    add_atom(sizes.xmax, sizes.ymax, zmin_down, 0);
}

// Determine whether an atom is near the edge of simulation box
const bool Bulk::on_edge(const double x, const double x_boundary) {
    return fabs(x - x_boundary) <= crys_struct.latconst / 2.1;
}

// =================================================
//            Implementation of Surface
// =================================================

// Constructor for Surface class
Surface::Surface(const double latconst, const int nnn) {
    init_statistics();
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

// Function to pick suitable extraction function
const void Surface::extract_surface(AtomReader* reader, const Femocs::SimuCell* cell) {
    if (reader->data.simu_type == "md")
        coordination_extract(reader, cell);
    else if (reader->data.simu_type == "kmc") kmc_extract(reader, cell);
}

// Extract surface by coordination analysis
const void Surface::kmc_extract(AtomReader* reader, const Femocs::SimuCell* cell) {
    int i;
    int N = reader->get_n_atoms();
    vector<bool> is_surface(N);

    // Get number and locations of surface atoms
    for (i = 0; i < N; ++i)
        is_surface[i] = (reader->get_type(i) == cell->type_surf);

    // Preallocate memory for Surface atoms
    reserve(accumulate(is_surface.begin(), is_surface.end(), 0));

    // Add surface atoms to Surface
    for (i = 0; i < N; ++i)
        if (is_surface[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i), reader->get_coord(i));

    init_statistics();
    calc_statistics();
}

// Extract surface by coordination analysis
const void Surface::coordination_extract(AtomReader* reader, const Femocs::SimuCell* cell) {
    int N = reader->get_n_atoms();
    int i, j;
    double r2, dx, dy, dz;
    double zmin = cell->zmin + crys_struct.latconst;

    vector<bool> is_surface(N);

    for (i = 0; i < N; ++i)
        is_surface[i] = (reader->get_z(i) > zmin) && (reader->get_coord(i) > 0)
                && (reader->get_coord(i) < crys_struct.nnn);

    // Preallocate memory for Surface atoms
    reserve(accumulate(is_surface.begin(), is_surface.end(), 0));

    for (i = 0; i < N; ++i)
        if (is_surface[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i), reader->get_coord(i));

    init_statistics();
    calc_statistics();
}

// =================================================
//              Implementation of Edge
// =================================================

// Constructor for Edge class
Edge::Edge(const double latconst, const int nnn) {
    init_statistics();
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

// Determine whether an atom is near the edge of simulation box
const bool Edge::on_edge(const double x, const double x_boundary) {
    return fabs(x - x_boundary) <= crys_struct.latconst / 2.1;
}

// Exctract the atoms near the simulation box sides
const void Edge::extract_edge(Surface* atoms, const Femocs::SimuCell* cell) {
    int i;
    int N = atoms->get_n_atoms();
    this->reserve(N);

    for (i = 0; i < N; ++i) {
        if (on_edge(atoms->get_x(i), cell->xmin))
            add_atom(cell->xmin, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_x(i), cell->xmax))
            add_atom(cell->xmax, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_y(i), cell->ymin))
            add_atom(atoms->get_x(i), cell->ymin, atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_y(i), cell->ymax))
            add_atom(atoms->get_x(i), cell->ymax, atoms->get_z(i), atoms->get_coordination(i));
    }

    init_statistics();
    calc_statistics();
}

} /* namespace femocs */
