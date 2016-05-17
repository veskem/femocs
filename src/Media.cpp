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
const void Bulk::rectangularize(const AtomReader::Sizes* sizes, const double r_cut) {
    double zmin_up = this->sizes.zmin + crys_struct.latconst / 2.1;
    double zmin_down = this->sizes.zmin;

    for (int i = 0; i < get_n_atoms(); ++i) {
        // Flatten the atoms on bottom layer
        if ((get_z(i) >= zmin_down) && (get_z(i) <= zmin_up)) set_z(i, zmin_down);

        Point2d point = get_point2d(i);
        int near1 = point.near(sizes->xmin, sizes->ymin, r_cut);
        int near2 = point.near(sizes->xmax, sizes->ymax, r_cut);

        // Flatten the atoms on the sides of simulation box
        if (near1 == 0) set_x(i, sizes->xmin);
        else if (near1 == 1) set_y(i, sizes->ymin);
        else if (near2 == 0) set_x(i, sizes->xmax);
        else if (near2 == 1) set_y(i, sizes->ymax);
    }
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
    double cutoff;
    int i, j;

    int n_atoms = get_n_atoms();

    vector<bool> is_dense(n_atoms);
    vector<int> coarse_indxs;

    Point2d origin2d((sizes.xmax + sizes.xmin) / 2, (sizes.ymax + sizes.ymin) / 2);

    // Create a map from atoms in- and outside the dense region
    for (i = 0; i < n_atoms; ++i)
        is_dense[i] = ((get_point2d(i).distance(origin2d) < r_cut) || (get_z(i) >= zmax));

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

        //double cutoff = 0.5*crys_struct.latconst * log( 10 + (point1.distance(origin3d) - r_cut) );
        cutoff = 0.3 * crys_struct.latconst * sqrt(1 + (point1.distance(origin3d) - r_cut));

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
    coarse_surf.reserve(n_coarse_atoms);
    for (int ci : coarse_indxs)
        if (ci >= 0) coarse_surf.add_atom(x[ci], y[ci], z[ci], coordination[ci]);

    coarse_surf.rectangularize(ar_sizes, 0.5 * crys_struct.latconst);

    Edge edge(crys_struct.latconst, crys_struct.nnn);
    edge.extract_edge(this, ar_sizes, cutoff);
    edge = edge.clean();
    edge = edge.coarsen(cutoff);
    //Edge edge2 = edge.coarsen(cutoff);

    dense_surf += coarse_surf;
    dense_surf += edge;

    return dense_surf;
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
        Point2d point = get_point2d(i);
        int near1 = point.near(sizes->xmin, sizes->ymin, r_cut);
        int near2 = point.near(sizes->xmax, sizes->ymax, r_cut);

        // Flatten the atoms on the sides of simulation box
        if (near1 == 0) set_x(i, sizes->xmin);
        else if (near1 == 1) set_y(i, sizes->ymin);
        else if (near2 == 0) set_x(i, sizes->xmax);
        else if (near2 == 1) set_y(i, sizes->ymax);
    }
}

const Surface Surface::rectangularize_vol2(const AtomReader::Sizes* sizes, const double r_cut) {
    // Loop through all the atoms
    for (int i = 0; i < get_n_atoms(); ++i) {
        Point2d point = get_point2d(i);
        int near1 = point.near(sizes->xmin, sizes->ymin, r_cut);
        int near2 = point.near(sizes->xmax, sizes->ymax, r_cut);

        // Flatten the atoms on the sides of simulation box
        if (near1 == 0) set_x(i, sizes->xmin);
        if (near1 == 1) set_y(i, sizes->ymin);
        if (near2 == 0) set_x(i, sizes->xmax);
        if (near2 == 1) set_y(i, sizes->ymax);
    }

    // Add atom to the lateral corners of surface
    add_atom(sizes->xmin, sizes->ymin, this->sizes.zmin, 0);
    add_atom(sizes->xmin, sizes->ymax, this->sizes.zmin, 0);
    add_atom(sizes->xmax, sizes->ymin, this->sizes.zmin, 0);
    add_atom(sizes->xmax, sizes->ymax, this->sizes.zmin, 0);

    // Remove the overlapping atoms
    return clean();
}

// Function to clean the Surface from overlapping atoms
const Surface Surface::clean() {
    Surface surf(crys_struct.latconst, crys_struct.nnn);
    int i, j;
    Point3d point1;
    int n_atoms = get_n_atoms();

    vector<bool> do_delete(n_atoms, false);

    // Loop through all the atoms
    for(i = 0; i < n_atoms-1; ++i) {
        // Skip already deleted atoms
        if (do_delete[i]) continue;

        point1 = get_point(i);
        for (j = i+1; j < n_atoms; ++j) {
            // Skip already deleted atoms
            if (do_delete[j]) continue;
            if (point1.distance2(get_point(j)) == 0)
                do_delete[j] = true;
        }
    }

    surf.reserve( n_atoms - accumulate(do_delete.begin(), do_delete.end(), 0) );
    for (i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            surf.add_atom(x[i], y[i], z[i], coordination[i]);

    surf.calc_statistics();
    return surf;
}

const bool Surface::contains(Point2d &p) {
    for (int i = 0; i < get_n_atoms(); ++i)
        if ( p.distance2(get_point2d(i)) == 0 )
            return true;

    return false;
}

const bool Surface::contains(Point3d &p) {
    for (int i = 0; i < get_n_atoms(); ++i)
        if ( p.distance(get_point(i)) == 0 )
            return true;
    return false;
}

// =================================================
//              Implementation of Edge
// =================================================

// Constructors for Edge class
Edge::Edge(const double latconst, const int nnn) :
        Medium() {
    crys_struct.latconst = latconst;
    crys_struct.nnn = nnn;
}
;
Edge::Edge() : Medium() {
    crys_struct.latconst = 0;
    crys_struct.nnn = 0;
}
;

// Exctract the atoms near the simulation box sides
const void Edge::extract_edge(Medium* atoms, const AtomReader::Sizes* sizes, const double r_cut) {
    int i;
    int N = atoms->get_n_atoms();

    // Reserve memory for atoms
    reserve(4 + N);

    // Add 4 atoms to the bottom corners of the edge
    add_atom(sizes->xmin, sizes->ymin, atoms->sizes.zmin, 0);
    add_atom(sizes->xmin, sizes->ymax, atoms->sizes.zmin, 0);
    add_atom(sizes->xmax, sizes->ymin, atoms->sizes.zmin, 0);
    add_atom(sizes->xmax, sizes->ymax, atoms->sizes.zmin, 0);

    for (i = 0; i < N; ++i) {
        Point2d point = atoms->get_point2d(i);
        int near1 = point.near(sizes->xmin, sizes->ymin, r_cut);
        int near2 = point.near(sizes->xmax, sizes->ymax, r_cut);

        if (near1 == 0)
            add_atom(sizes->xmin, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        else if (near1 == 1)
            add_atom(atoms->get_x(i), sizes->ymin, atoms->get_z(i), atoms->get_coordination(i));
        else if (near2 == 0)
            add_atom(sizes->xmax, atoms->get_y(i), atoms->get_z(i), atoms->get_coordination(i));
        else if (near2 == 1)
            add_atom(atoms->get_x(i), sizes->ymax, atoms->get_z(i), atoms->get_coordination(i));
    }

    calc_statistics();
}

const Edge Edge::coarsen(const double r_cut) {
    int i, j;
    int n_atoms = get_n_atoms();

    Point3d point1;

    Edge edge(crys_struct.latconst, crys_struct.nnn);
    vector<bool> is_dense(n_atoms, false);

    // Mark atoms that are too close and need to be deleted
    for (i = 0; i < (n_atoms - 1); ++i) {
        // Skip already deleted atoms
        if (is_dense[i]) continue;

        // Mark too close nodes
        point1 = get_point(i);
        for (j = i + 1; j < n_atoms; ++j) {
            // Skip already deleted atoms
            if (is_dense[j]) continue;

            if (point1.distance(get_point(j)) < r_cut) is_dense[j] = true;
        }
    }

    // Reserve memory for edge atoms
    edge.reserve(n_atoms - accumulate(is_dense.begin(), is_dense.end(), 0));
    // Compile coarsened Edge
    for (i = 0; i < n_atoms; ++i)
        if (!is_dense[i]) edge.add_atom(get_x(i), get_y(i), get_z(i), get_coordination(i));

    edge.calc_statistics();
    return edge;

}

const Edge Edge::clean() {
    Edge edge(crys_struct.latconst, crys_struct.nnn);
    int i, j;
    Point3d point1;
    int n_atoms = get_n_atoms();

    vector<bool> do_delete(n_atoms, false);

    // Loop through all the atoms
    for(i = 0; i < n_atoms-1; ++i) {
        // Skip already deleted atoms
        if (do_delete[i]) continue;

        point1 = get_point(i);
        for (j = i+1; j < n_atoms; ++j) {
            // Skip already deleted atoms
            if (do_delete[j]) continue;
            if (point1.distance2(get_point(j)) == 0)
                do_delete[j] = true;
        }
    }

    edge.reserve( n_atoms - accumulate(do_delete.begin(), do_delete.end(), 0) );
    for (i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            edge.add_atom(x[i], y[i], z[i], coordination[i]);

    edge.calc_statistics();
    return edge;
}

} /* namespace femocs */
