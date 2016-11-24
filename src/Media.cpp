/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Media.h"
#include <numeric>

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
    add_atom( Atom(-1, Point3(sizes->xmax, sizes->ymax, sizes->zmaxbox), 0) );
    add_atom( Atom(-1, Point3(sizes->xmax, sizes->ymin, sizes->zmaxbox), 0) );
    add_atom( Atom(-1, Point3(sizes->xmin, sizes->ymax, sizes->zmaxbox), 0) );
    add_atom( Atom(-1, Point3(sizes->xmin, sizes->ymin, sizes->zmaxbox), 0) );

    init_statistics();
    calc_statistics();
}

// =================================================
//              Implementation of Bulk
// =================================================

// Constructor for Bulk
Bulk::Bulk() : Medium() {};

// Function to make bulk material with nodes on surface and on by-hand made bottom coordinates
const void Bulk::generate_simple(const AtomReader::Sizes* sizes) {
    int M = 4; // total number of nodes
    reserve(M);

    // Add atoms to the bottom corner edges of simulation cell
    add_atom( Atom(-1, Point3(sizes->xmin, sizes->ymin, sizes->zminbox), TYPES.BULK) );
    add_atom( Atom(-1, Point3(sizes->xmin, sizes->ymax, sizes->zminbox), TYPES.BULK) );
    add_atom( Atom(-1, Point3(sizes->xmax, sizes->ymin, sizes->zminbox), TYPES.BULK) );
    add_atom( Atom(-1, Point3(sizes->xmax, sizes->ymax, sizes->zminbox), TYPES.BULK) );

    init_statistics();
    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::extract(AtomReader* reader) {
    const int N = reader->get_n_atoms();
    vector<bool> is_bulk(N);

    for (int i = 0; i < N; ++i)
        is_bulk[i] = reader->get_type(i) != TYPES.VACANCY;

    reserve(vector_sum(is_bulk));

    for (int i = 0; i < N; ++i)
        if (is_bulk[i])
            add_atom(reader->get_atom(i));

    calc_statistics();
}

// Function to extract bulk material from input atomistic data
const void Bulk::rectangularize(const AtomReader::Sizes* sizes, const double eps, const double latconst) {
    const double zmin_up = this->sizes.zmin + latconst / 2.1;
    const double zmin_down = this->sizes.zmin;

    for (int i = 0; i < get_n_atoms(); ++i) {
        Point3 point = get_point(i);

        // Flatten the atoms on bottom layer
        if ( point.z <= zmin_up ) set_z(i, zmin_down);

        // Flatten the atoms on the sides of simulation box
        if ( on_boundary(point.x, sizes->xmin, eps) ) set_x(i, sizes->xmin);
        else if ( on_boundary(point.y, sizes->ymin, eps) ) set_y(i, sizes->ymin);
        else if ( on_boundary(point.x, sizes->xmax, eps) ) set_x(i, sizes->xmax);
        else if ( on_boundary(point.y, sizes->ymax, eps) ) set_y(i, sizes->ymax);
    }
}

// =================================================
//            Implementation of Surface
// =================================================

// Constructors for Surface class
Surface::Surface() : Medium() {}

Surface::Surface(const int n_atoms) : Medium(n_atoms) {}

// Generate Surface with 4 atoms at the corners and 4 at the middle edges of simulation cell
const void Surface::generate_simple(const AtomReader::Sizes* sizes, const double z) {
    // Reserve memory for atoms
    reserve(4);

    // Add 4 atoms to the corners of simulation cell
    // The insertion order determines the orientation of big surface triangles
    add_atom( Point3(sizes->xmin, sizes->ymin, z) );
    add_atom( Point3(sizes->xmax, sizes->ymin, z) );
    add_atom( Point3(sizes->xmax, sizes->ymax, z) );
    add_atom( Point3(sizes->xmin, sizes->ymax, z) );
}

// Generate edge with regular atom distribution between surface corners
const void Surface::generate_middle(const AtomReader::Sizes* sizes, const double z, const double r_cut) {
    const int n_atoms_per_side_x = sizes->xbox / r_cut + 1;
    const int n_atoms_per_side_y = sizes->ybox / r_cut + 1;

    // Reserve memory for atoms
    reserve( 2 * (n_atoms_per_side_x + n_atoms_per_side_y) - 4 );

    // Add atoms in x-direction
    for (int i = 1; i < n_atoms_per_side_x; ++i) {
        double coordinate = sizes->xmin + i * sizes->xbox / n_atoms_per_side_x;
        add_atom( Atom(-1, Point3(sizes->xmin, coordinate, z), 0) );
        add_atom( Atom(-1, Point3(sizes->xmax, coordinate, z), 0) );
    }

    // Add atoms in y-direction
    for (int i = 1; i < n_atoms_per_side_y; ++i) {
        double coordinate = sizes->ymin + i * sizes->ybox / n_atoms_per_side_y;
        add_atom( Atom(-1, Point3(coordinate, sizes->ymin, z), 0) );
        add_atom( Atom(-1, Point3(coordinate, sizes->ymax, z), 0) );
    }
}

// Extract surface by the atom types
const void Surface::extract(AtomReader* reader) {
    const int n_atoms = reader->get_n_atoms();
    vector<bool> is_surface(n_atoms);

    // Get number and locations of surface atoms
    for (int i = 0; i < n_atoms; ++i)
        is_surface[i] = reader->get_type(i) == TYPES.SURFACE;

    // Preallocate memory for Surface atoms
    reserve(vector_sum(is_surface));

    // Add surface atoms to Surface
    for (int i = 0; i < n_atoms; ++i)
        if (is_surface[i])
            add_atom(reader->get_atom(i));
}

// Function to coarsen the atoms with coarsener
const Surface Surface::coarsen(Coarseners &coarseners, const AtomReader::Sizes* reader) {
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);

    Surface corners; corners.generate_simple(reader, coarseners.zmean);
    Surface middle; middle.generate_middle(reader, coarseners.zmean, coarseners.r0_inf);
    middle.sort_atoms(3, "down", origin2d);
    this->sort_atoms(3, 2, "down", origin2d);

    Surface union_surf;
    union_surf += corners;
    union_surf += middle;
    union_surf.add(this);

    return union_surf.clean(coarseners);
}

// Function to flatten the atoms on the sides of simulation box
const Surface Surface::rectangularize(const AtomReader::Sizes* sizes, const double eps, const double latconst) {
    Coarseners coarseners;
    coarseners.attach_coarsener(make_shared<ConstCoarsener>(latconst / 3.0));

    Edge edge;
    edge.extract(this, sizes, eps);

    Surface surf;
    surf += edge;
    surf.add(this);

    return surf.clean(coarseners);
}

// Clean the surface from atoms that are too close to each other
const Surface Surface::clean(Coarseners &coarseners) {
    const int n_atoms = get_n_atoms();
    vector<bool> do_delete(n_atoms, false);

    // Loop through all the atoms
    for(int i = 0; i < n_atoms-1; ++i) {
        // Skip already deleted atoms
        if (do_delete[i]) continue;

        Point3 point1 = get_point(i);
        coarseners.pick_cutoff(point1);

        for (int j = i+1; j < n_atoms; ++j) {
            // Skip already deleted atoms
            if (do_delete[j]) continue;
            do_delete[j] = coarseners.nearby(point1, get_point(j));
        }
    }

    Surface surf( n_atoms - vector_sum(do_delete) );
    for (int i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            surf.add_atom(get_atom(i));

    surf.calc_statistics();

    return surf;
}

// Function to delete atoms that are separate from others
// Atom is considered lonely if its coordination is lower than coord_min
const Surface Surface::clean_lonely_atoms(const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const int n_atoms = get_n_atoms();
    const int coord_min = 2;

    vector<bool> atom_not_lonely(n_atoms, false);
    vector<int> coords(n_atoms, 0);

    // Loop through all the atoms
    for(int i = 0; i < n_atoms-1; ++i) {
        // Skip already processed atoms
        if (atom_not_lonely[i]) continue;

        Point3 point1 = get_point(i);

        // Loop through the possible neighbours of the atom
        for (int j = i+1; j < n_atoms; ++j) {
            // Skip already processed atoms
            if (atom_not_lonely[j]) continue;

            if (point1.distance2(get_point(j)) <= r_cut2) {
                atom_not_lonely[i] = ++coords[i] >= coord_min;
                atom_not_lonely[j] = ++coords[j] >= coord_min;
            }
        }
    }

    Surface surf( vector_sum(atom_not_lonely) );
    for (int i = 0; i < n_atoms; ++i)
        if (atom_not_lonely[i])
            surf.add_atom(get_atom(i));

    surf.calc_statistics();
    return surf;
}

inline double Surface::smooth_function(double distance, double smooth_factor) const {
    const double a = 1.0;
    return a * exp(-1.0 * distance / smooth_factor);
}

const void Surface::smoothen(double radius, double smooth_factor, double r_cut) {
    // Calculate the horizontal span of the surface
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);
    smoothen(origin2d, radius, smooth_factor, r_cut);
}

const void Surface::smoothen(const Point2 &origin, double radius, double smooth_factor, double r_cut) {
    if (smooth_factor < 0.01) return;

    const int n_atoms = get_n_atoms();
    const double radius2 = radius * radius;

    // Find atoms in interesting region
    vector<bool> hot_atom(n_atoms);
    for(int i = 0; i < n_atoms; ++i)
        hot_atom[i] = origin.distance2(get_point2(i)) <= radius2;

    const int n_smooth = accumulate(hot_atom.begin(), hot_atom.end(), 0);

    // Transfer the points in interesting region into new surface
    Surface temp_surf(n_smooth);
    for(int i = 0; i < n_atoms; ++i)
        if(hot_atom[i])
            temp_surf.add_atom(get_point(i));

    // Smoothen the surface
    temp_surf.smoothen(smooth_factor, r_cut);

    // Write smoothened points back to surface
    int j = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (hot_atom[i])
            set_point(i, temp_surf.get_point(j++));
}

const void Surface::smoothen(double smooth_factor, double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const int n_atoms = get_n_atoms();

    // Make copy of points so that the old positions won't interfere with already smoothed ones
    vector<Point3> points; points.reserve(n_atoms);
    for(int i = 0; i < n_atoms; ++i)
        points.push_back(get_point(i));

    // Vector for sum of weights
    vector<double> weights_sum(n_atoms, 1.0);

    // Smooth the vertices
    for (int i = 0; i < n_atoms-1; ++i) {
        Point3 point1 = points[i];

        for (int j = i+1; j < n_atoms; ++j) {
            Point3 point2 = points[j];
            double distance2 = point1.distance2(point2);
            if (distance2 > r_cut2) continue;

            double weight = smooth_function(sqrt(distance2), smooth_factor);
            atoms[i].point += point2 * weight;
            atoms[j].point += point1 * weight;
            weights_sum[i] += weight;
            weights_sum[j] += weight;
        }
    }

    // Normalise smoothed vertices
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].point *= 1.0 / weights_sum[i];
}

// =================================================
//              Implementation of Edge
// =================================================

// Constructors for Edge class
Edge::Edge() : Medium() {};

// Exctract the atoms near the simulation box sides
const void Edge::extract(const Medium* atoms, const AtomReader::Sizes* sizes, const double eps) {
    const int n_atoms = atoms->get_n_atoms();

    // Reserve memory for atoms
    reserve(n_atoms);

    // Get the atoms from edge areas
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms->get_point(i);
        const bool near1 = on_boundary(point.x, sizes->xmin, sizes->xmax, eps);
        const bool near2 = on_boundary(point.y, sizes->ymin, sizes->ymax, eps);

        if (near1 || near2)
            add_atom(atoms->get_atom(i));
    }

    // Loop through all the added atoms and flatten the atoms on the sides of simulation box
    for (int i = 0; i < get_n_atoms(); ++i) {
        Point3 point = get_point(i);
        if ( on_boundary(point.x, sizes->xmin, eps) ) set_x(i, sizes->xmin);
        if ( on_boundary(point.x, sizes->xmax, eps) ) set_x(i, sizes->xmax);
        if ( on_boundary(point.y, sizes->ymin, eps) ) set_y(i, sizes->ymin);
        if ( on_boundary(point.y, sizes->ymax, eps) ) set_y(i, sizes->ymax);
    }
}

} /* namespace femocs */
