/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Media.h"
#include <algorithm>

using namespace std;
namespace femocs {

// =================================================
//             Implementation of Vacuum
// =================================================

// Constructor of Vacuum class
Vacuum::Vacuum() :
        Medium() {
}

// Generates Vacuum by adding four points to the top of simulation cell
const void Vacuum::generate_simple(const AtomReader::Sizes* sizes) {
    int M = 4; // total number of nodes
    reserve(M);

    // Add points to the xy-plane corners on current layer
    add_atom(sizes->xmax, sizes->ymax, sizes->zmaxbox, 0);
    add_atom(sizes->xmax, sizes->ymin, sizes->zmaxbox, 0);
    add_atom(sizes->xmin, sizes->ymax, sizes->zmaxbox, 0);
    add_atom(sizes->xmin, sizes->ymin, sizes->zmaxbox, 0);

    init_statistics();
    calc_statistics();
}

// =================================================
//              Implementation of Bulk
// =================================================

// Constructor for Bulk class
Bulk::Bulk(const double latconst, const int nnn) :
        Medium() {
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const void Bulk::extract_reduced_bulk(Surface* surf, const AtomReader::Sizes* sizes) {
    int i;
//    shared_ptr<Edge> edge(new Edge(crys_struct.latconst, crys_struct.nnn));
//    edge->extract_edge(surf, sizes);

    int n_surf = surf->get_n_atoms();
    int n_edge = 4; //edge->getN();

    // Reserve memory for bulk atoms and their parameters
    reserve(n_surf + n_edge);

    for (i = 0; i < n_surf; ++i)
        add_atom(surf->get_x(i), surf->get_y(i), surf->get_z(i), surf->get_coordination(i));

//    for (i = 0; i < n_edge; ++i)
//        add_atom(edge->get_x(i), edge->get_y(i), sizes->zmin, edge->get_coordination(i));

    // Add extra atoms to the bottom corner edges of simulation cell
    add_atom(sizes->xmin, sizes->ymin, sizes->zminbox, 0);
    add_atom(sizes->xmin, sizes->ymax, sizes->zminbox, 0);
    add_atom(sizes->xmax, sizes->ymin, sizes->zminbox, 0);
    add_atom(sizes->xmax, sizes->ymax, sizes->zminbox, 0);

    calc_statistics();
}

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const void Bulk::generate_simple(const AtomReader::Sizes* sizes, const AtomReader::Types* types) {
    int M = 4; // total number of nodes
    reserve(M);

    // Add atoms to the bottom corner edges of simulation cell
    add_atom(sizes->xmin, sizes->ymin, sizes->zminbox, types->type_bulk);
    add_atom(sizes->xmin, sizes->ymax, sizes->zminbox, types->type_bulk);
    add_atom(sizes->xmax, sizes->ymin, sizes->zminbox, types->type_bulk);
    add_atom(sizes->xmax, sizes->ymax, sizes->zminbox, types->type_bulk);

    init_statistics();
    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract_truncated_bulk_old(AtomReader* reader) {
    int i;
    const int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);
    vector<bool> is_surf(N);

    // Get number and locations of bulk atoms
    for (i = 0; i < N; ++i)
        is_bulk[i] = (reader->get_z(i) > reader->sizes.zminbox)
                && (reader->get_type(i) == reader->types.type_bulk
                        || reader->get_type(i) == reader->types.type_surf);

    for (i = 0; i < N; ++i)
        is_surf[i] = (reader->get_coordination(i) > 0)
                && (reader->get_coordination(i) < crys_struct.nnn);

    reserve(accumulate(is_bulk.begin(), is_bulk.end(), 0));

    // Add bulk atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && !is_surf[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));

    // Add surface atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i] && is_surf[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract_truncated_bulk(AtomReader* reader) {
    int i;
    const int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);

    // Get number and locations of bulk atoms
    for (i = 0; i < N; ++i)
        is_bulk[i] = (reader->get_z(i) > reader->sizes.zminbox)
                && (reader->get_type(i) == reader->types.type_bulk);

    reserve(accumulate(is_bulk.begin(), is_bulk.end(), 0));

    // Add bulk atoms
    for (i = 0; i < N; ++i)
        if (is_bulk[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract_bulk(AtomReader* reader) {
    int i;
    int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);

    for (i = 0; i < N; ++i)
        is_bulk[i] = reader->get_type(i) != reader->types.type_vacancy;

    reserve(accumulate(is_bulk.begin(), is_bulk.end(), 0));

    for (i = 0; i < N; ++i)
        if (is_bulk[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::rectangularize(const AtomReader::Sizes* sizes) {
    double zmin_up = this->sizes.zmin + crys_struct.latconst / 2.1;
    double zmin_down = this->sizes.zmin;

    for (int i = 0; i < get_n_atoms(); ++i) {
        // Flatten the atoms on bottom layer
        if ((get_z(i) >= zmin_down) && (get_z(i) <= zmin_up)) set_z(i, zmin_down);

        // Flatten the atoms on the sides of simulation box
        if (on_edge(get_x(i), sizes->xmin)) set_x(i, sizes->xmin);
        if (on_edge(get_x(i), sizes->xmax)) set_x(i, sizes->xmax);
        if (on_edge(get_y(i), sizes->ymin)) set_y(i, sizes->ymin);
        if (on_edge(get_y(i), sizes->ymax)) set_y(i, sizes->ymax);
    }

    // Add atoms to the bottom corner of the simulation cell
//    add_atom(this->sizes.xmin, this->sizes.ymin, zmin_down, 0);
//    add_atom(this->sizes.xmin, this->sizes.ymax, zmin_down, 0);
//    add_atom(this->sizes.xmax, this->sizes.ymin, zmin_down, 0);
//    add_atom(this->sizes.xmax, this->sizes.ymax, zmin_down, 0);
}

// Determine whether an atom is near the edge of simulation box
const bool Bulk::on_edge(const double x, const double x_boundary) {
    return fabs(x - x_boundary) <= crys_struct.latconst / 2.0;
}

// =================================================
//            Implementation of Surface
// =================================================

// Constructor for Surface class
Surface::Surface(const double latconst, const int nnn) :
        Medium() {
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

const Surface Surface::coarsen(const double r_cut, const double zmax, const AtomReader::Sizes* ar_sizes) {
    Surface dense_surf(crys_struct.latconst, crys_struct.nnn);
    Surface coarse_surf(crys_struct.latconst, crys_struct.nnn);
    const double r_cut2 = r_cut * r_cut;
    int i, j;

    int n_atoms = get_n_atoms();

    vector<bool> is_dense(n_atoms);
    vector<int> coarse_indxs;

    Point2d origin2d((sizes.xmax + sizes.xmin) / 2, (sizes.ymax + sizes.ymin) / 2);

    // Create a map from atoms in- and outside the dense region
    for (i = 0; i < n_atoms; ++i)
        is_dense[i] = ( (get_point2d(i).distance(origin2d) < r_cut) || (get_z(i) >= zmax) );

    int n_dense_atoms = accumulate(is_dense.begin(), is_dense.end(), 0);
    int n_coarse_atoms = n_atoms - n_dense_atoms;

    dense_surf.reserve(n_dense_atoms);
    coarse_indxs.reserve(n_coarse_atoms);

    // Copy the all the atoms from dense region and save the indices of atoms to be coarsened
    for (i = 0; i < n_atoms; ++i) {
        if (is_dense[i])
            dense_surf.add_atom(x[i], y[i], z[i], coordination[i]);
        else
            coarse_indxs.push_back(i);
    }

    Point3d origin3d(origin2d[0], origin2d[1], zmax);

    // Mark the atoms that are too close to each other in coarse region
    for (i = 0; i < (coarse_indxs.size() - 1); ++i) {
        // Skip already deleted atoms
        if (coarse_indxs[i] < 0) continue;

        Point3d point1 = get_point(coarse_indxs[i]);

        double dist = point1.distance(origin3d);
        //double cutoff = 0.5*crys_struct.latconst * log( 10 + (dist - r_cut) );
        double cutoff = 0.5*crys_struct.latconst * sqrt( 1 + (dist - r_cut) );

        for (j = i + 1; j < coarse_indxs.size(); ++j) {
            // Skip already deleted atoms
            if (coarse_indxs[j] < 0) continue;

            // Delete atom that is too close to given atom
            if (point1.distance(get_point(coarse_indxs[j])) < cutoff) {
                coarse_indxs[j] = -1;
                --n_coarse_atoms;
            }
        }
    }

    require(n_coarse_atoms >= 0, "Invalid number of coarse atoms!");

    // Compile the sub-surface from coarsened atoms
    coarse_surf.reserve(4 + n_coarse_atoms);

    for (int ci : coarse_indxs)
        if (ci >= 0) coarse_surf.add_atom(x[ci], y[ci], z[ci], coordination[ci]);

    coarse_surf.rectangularize(ar_sizes, 1.0*crys_struct.latconst);

    // Add 4 atoms to the corners of input surface
    coarse_surf.add_atom(sizes.xmin, sizes.ymin, sizes.zmin, 0);
    coarse_surf.add_atom(sizes.xmin, sizes.ymax, sizes.zmin, 0);
    coarse_surf.add_atom(sizes.xmax, sizes.ymin, sizes.zmin, 0);
    coarse_surf.add_atom(sizes.xmax, sizes.ymax, sizes.zmin, 0);

    return coarse_surf + dense_surf;
}


// Function to pick suitable extraction function
const void Surface::extract_surface(AtomReader* reader) {
    //if (reader->types.simu_type == "md")
    //    extract_by_coordination(reader);
    //else if (reader->types.simu_type == "kmc")
    extract_by_type(reader);

    calc_statistics();
}

// Extract surface by coordination analysis
const void Surface::extract_by_type(AtomReader* reader) {
    int i;
    int N = reader->get_n_atoms();
    vector<bool> is_surface(N);

    // Get number and locations of surface atoms
    for (i = 0; i < N; ++i)
        is_surface[i] = (reader->get_type(i) == reader->types.type_surf);

    // Preallocate memory for Surface atoms
    reserve(accumulate(is_surface.begin(), is_surface.end(), 0));

    // Add surface atoms to Surface
    for (i = 0; i < N; ++i)
        if (is_surface[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));
}

// Extract surface by coordination analysis
const void Surface::extract_by_coordination(AtomReader* reader) {
    int n_atoms = reader->get_n_atoms();
    int i;
    double zmin = reader->sizes.zminbox + crys_struct.latconst;

    vector<bool> is_surface(n_atoms);

    for (i = 0; i < n_atoms; ++i)
        is_surface[i] = (reader->get_z(i) > zmin) && (reader->get_coordination(i) > 0)
                && (reader->get_coordination(i) < crys_struct.nnn);

    // Preallocate memory for Surface atoms
    reserve(accumulate(is_surface.begin(), is_surface.end(), 0));

    for (i = 0; i < n_atoms; ++i)
        if (is_surface[i])
            add_atom(reader->get_x(i), reader->get_y(i), reader->get_z(i),
                    reader->get_coordination(i));
}

// Function to flatten the atoms on the sides of simulation box
const void Surface::rectangularize(const AtomReader::Sizes* sizes, const double r_cut) {
    // Loop through all the atoms
    for (int i = 0; i < get_n_atoms(); ++i) {
        // Flatten the atoms on the sides of simulation box
        if (on_edge(get_x(i), sizes->xmin, r_cut)) set_x(i, sizes->xmin);
        if (on_edge(get_x(i), sizes->xmax, r_cut)) set_x(i, sizes->xmax);
        if (on_edge(get_y(i), sizes->ymin, r_cut)) set_y(i, sizes->ymin);
        if (on_edge(get_y(i), sizes->ymax, r_cut)) set_y(i, sizes->ymax);
    }
}

// Determine whether an atom is near the edge of simulation box
const bool Surface::on_edge(const double x, const double x_boundary, const double r_cut) {
    return fabs(x - x_boundary) <= r_cut;
}

// =================================================
//              Implementation of Edge
// =================================================

// Constructor for Edge class
Edge::Edge(const double latconst, const int nnn) :
        Medium() {
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;

// Determine whether an atom is near the edge of simulation box
const bool Edge::on_edge(const double x, const double x_boundary) {
    return fabs(x - x_boundary) <= crys_struct.latconst / 2.1;
}

// Exctract the atoms near the simulation box sides
const void Edge::extract_edge(Surface* atoms, const AtomReader::Sizes* sizes) {
    int i;
    int N = atoms->get_n_atoms();
    this->reserve(N);

    for (i = 0; i < N; ++i) {
        if (on_edge(atoms->get_x(i), sizes->xmin))
            add_atom(sizes->xmin, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_x(i), sizes->xmax))
            add_atom(sizes->xmax, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_y(i), sizes->ymin))
            add_atom(atoms->get_x(i), sizes->ymin, atoms->get_z(i), atoms->get_coordination(i));
        if (on_edge(atoms->get_y(i), sizes->ymax))
            add_atom(atoms->get_x(i), sizes->ymax, atoms->get_z(i), atoms->get_coordination(i));
    }

    calc_statistics();
}

} /* namespace femocs */
