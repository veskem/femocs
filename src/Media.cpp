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

Media::Media() : Medium() {}

Media::Media(const int n_atoms) : Medium(n_atoms) {}

Media::Media(const Medium::Sizes& ar_sizes, const double z) {
    generate_simple(ar_sizes, z);
}

// Generate Surface with 4 atoms at the corners and 4 at the middle edges of simulation cell
void Media::generate_simple(const Medium::Sizes& ar_sizes, const double z) {
    // Reserve memory for atoms
    reserve(4);

    // Add 4 atoms to the corners of simulation cell
    // The insertion order determines the orientation of big surface triangles
    append( Point3(ar_sizes.xmin, ar_sizes.ymin, z) );
    append( Point3(ar_sizes.xmax, ar_sizes.ymin, z) );
    append( Point3(ar_sizes.xmax, ar_sizes.ymax, z) );
    append( Point3(ar_sizes.xmin, ar_sizes.ymax, z) );

    calc_statistics();
}

// Generate edge with regular atom distribution between surface corners
void Media::generate_middle(const Medium::Sizes& s, const double z, const double dist) {
    require(dist > 0, "Invalid dist between atoms: " + to_string(dist));
    const int n_atoms_per_side_x = s.xbox / dist + 1;
    const int n_atoms_per_side_y = s.ybox / dist + 1;

    // Reserve memory for atoms
    reserve( 2 * (n_atoms_per_side_x + n_atoms_per_side_y) - 4 );

    // Add atoms on x-edge
    for (int i = 1; i < n_atoms_per_side_x; ++i) {
        double y = s.xmin + i * s.xbox / n_atoms_per_side_x;
        append( Point3(s.xmin, y, z) );
        append( Point3(s.xmax, y, z) );
    }

    // Add atoms on y-edge
    for (int i = 1; i < n_atoms_per_side_y; ++i) {
        double x = s.ymin + i * s.ybox / n_atoms_per_side_y;
        append( Point3(x, s.ymin, z) );
        append( Point3(x, s.ymax, z) );
    }
}

// Extract surface by the atom types
void Media::extract(const AtomReader& reader, const int type, const bool invert) {
    const int coord_min = 2;
    const int n_atoms = reader.size();
    vector<bool> is_type(n_atoms);

    // Get number and locations of atoms of desired type
    if (!invert) {
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = reader.get_marker(i) == type;
    } else {
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = reader.get_marker(i) != type;
    }

    // Clean lonely atoms; atom is considered lonely if its coordination is lower than coord_min
    if (reader.get_nborlist_size() == n_atoms)
        for (int i = 0; i < n_atoms; ++i)
            if (is_type[i]) {
                int n_nbors = 0;
                for (int nbor : reader.get_neighbours(i)) {
                    require(nbor >= 0 && nbor < n_atoms, "Invalid index: " + to_string(nbor));
                    if (is_type[nbor]) n_nbors++;
                }

                is_type[i] = n_nbors >= coord_min;
            }

    // Preallocate memory for atoms
    reserve(vector_sum(is_type));

    // Store the atoms
    for (int i = 0; i < n_atoms; ++i)
        if (is_type[i])
            append(reader.get_atom(i));
            
    calc_statistics();        
}

// Function extend the flat area by generating additional atoms
Media Media::extend(const string &file_name, Coarseners &coarseners) {
    AtomReader reader;
    reader.import_file(file_name);

    Media stretched(reader.size());
    stretched += reader;
    stretched.calc_statistics();
    stretched.sizes.zmean = stretched.sizes.zmin;

    return stretched.coarsen(coarseners);
}

// Function extend the flat area by generating additional atoms
Media Media::extend(const double latconst, const double box_width, Coarseners &coarseners) {
    const int n_atoms = size();

    calc_statistics();
    const double current_box_width = min(sizes.xbox, sizes.ybox);
    const double desired_box_width = box_width * sizes.zbox;

    // get over estimation for number of generated points
    const int n_generated = pow(desired_box_width / latconst + 1, 2) - pow(current_box_width / latconst - 1, 2);

    // copy input points without modification
    Media stretched( max(0, n_generated) );

    const double W = (desired_box_width - current_box_width) / 2.0;  // generation area width

    // if the input surface isn't sufficiently wide, add atoms to it, otherwise just add the boundary nodes
    if (W > 0) {
        // add points outside already existing area
        for (double y = sizes.ymin - W; y <= sizes.ymax + W; y += latconst)
            for (double x = sizes.xmin - W; x <= sizes.xmax + W; x += latconst)
                if ( x < sizes.xmin || x > sizes.xmax || y < sizes.ymin || y > sizes.ymax)
                    stretched.append( Point3(x, y, coarseners.centre.z) );
        stretched.calc_statistics();
    }
    else {
        stretched.copy_statistics(*this);
        stretched.sizes.zmean = coarseners.centre.z;
    }

    return stretched.coarsen(coarseners);
}

// Function to coarsen the atoms with coarsener
Media Media::coarsen(Coarseners &coarseners) {
    Media corners, middle, union_surf;

    corners.generate_simple(sizes, sizes.zmean);
    middle.generate_middle( sizes, sizes.zmean, coarseners.get_r0_inf(sizes) );
    sort_atoms(3, "down");

    union_surf += corners;
    union_surf += middle;
    union_surf.add(this);

    return union_surf.clean(coarseners);
}

// Clean the surface from atoms that are too close to each other
Media Media::get_nanotip(const double radius) {
    const int n_atoms = size();
    const double radius2 = radius*radius;

    vector<bool> do_delete;
    do_delete.reserve(n_atoms);

    Point2 centre(sizes.xmid, sizes.ymid);

    // Loop through all the atoms
    for (int i = 0; i < n_atoms; ++i)
       do_delete.push_back( centre.distance2(get_point2(i)) > radius2 );

    Media nanotip( n_atoms - vector_sum(do_delete) );
    for (int i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            nanotip.append(get_atom(i));

    nanotip.calc_statistics();
    return nanotip;
}

// Clean the surface from atoms that are too close to each other
Media Media::clean(Coarseners &coarseners) {
    const int n_atoms = size();
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

    Media surf( n_atoms - vector_sum(do_delete) );
    for (int i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            surf.append(get_atom(i));

    surf.calc_statistics();
    return surf;
}

// Remove the atoms that are too far from surface faces
void Media::clean(const TetgenMesh& mesh, const double r_cut) {
    if (r_cut <= 0) return;

    const int n_atoms = size();
    RaySurfaceIntersect rsi(&mesh);
    rsi.precompute_triangles();

    vector<Atom> _atoms;
    _atoms.reserve(n_atoms);

    for (int i = 0; i < n_atoms; ++i) {
        if ( rsi.near_surface(Vec3(get_point(i)), r_cut) )
            _atoms.push_back(atoms[i]);
    }

    atoms = _atoms;
    calc_statistics();
}

inline double Media::smooth_function(const double distance, const double smooth_factor) const {
    const double a = 1.0;
    return a * exp(-1.0 * distance / smooth_factor);
}

void Media::smoothen(const double radius, const double smooth_factor, const double r_cut) {
    // Calculate the horizontal span of the surface
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);
    smoothen(origin2d, radius, smooth_factor, r_cut);
}

void Media::smoothen(const Point2 &origin, const double radius, const double smooth_factor, const double r_cut) {
    if (smooth_factor < 0.01) return;

    const int n_atoms = size();
    const double radius2 = radius * radius;

    // Find atoms in interesting region
    vector<bool> hot_atom(n_atoms);
    for(int i = 0; i < n_atoms; ++i)
        hot_atom[i] = origin.distance2(get_point2(i)) <= radius2;

    const int n_smooth = accumulate(hot_atom.begin(), hot_atom.end(), 0);

    // Transfer the points in interesting region into new surface
    Media temp_surf(n_smooth);
    for(int i = 0; i < n_atoms; ++i)
        if(hot_atom[i])
            temp_surf.append(get_point(i));

    // Smoothen the surface
    temp_surf.smoothen(smooth_factor, r_cut);

    // Write smoothened points back to surface
    int j = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (hot_atom[i])
            set_point(i, temp_surf.get_point(j++));
}

void Media::smoothen(const double smooth_factor, const double r_cut) {
    const double r_cut2 = r_cut * r_cut;
    const int n_atoms = size();

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
void Edge::extract(const Medium* atoms, const AtomReader::Sizes& ar_sizes, const double eps) {
    const int n_atoms = atoms->size();

    // Reserve memory for atoms
    reserve(n_atoms);

    // Get the atoms from edge areas
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms->get_point(i);
        const bool near1 = on_boundary(point.x, ar_sizes.xmin, ar_sizes.xmax, eps);
        const bool near2 = on_boundary(point.y, ar_sizes.ymin, ar_sizes.ymax, eps);

        if (near1 || near2)
            append(atoms->get_atom(i));
    }

    // Loop through all the added atoms and flatten the atoms on the sides of simulation box
    for (int i = 0; i < size(); ++i) {
        Point3 point = get_point(i);
        if ( on_boundary(point.x, ar_sizes.xmin, eps) ) set_x(i, ar_sizes.xmin);
        if ( on_boundary(point.x, ar_sizes.xmax, eps) ) set_x(i, ar_sizes.xmax);
        if ( on_boundary(point.y, ar_sizes.ymin, eps) ) set_y(i, ar_sizes.ymin);
        if ( on_boundary(point.y, ar_sizes.ymax, eps) ) set_y(i, ar_sizes.ymax);
    }
}

} /* namespace femocs */
