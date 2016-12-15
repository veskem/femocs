/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Media.h"
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>  

using namespace std;
namespace femocs {

Media::Media() : Medium() {}

Media::Media(const int n_atoms) : Medium(n_atoms) {}


Media::Media(const Medium::Sizes& ar_sizes, const double z) {
    generate_simple(ar_sizes, z);
}

// Generate Surface with 4 atoms at the corners and 4 at the middle edges of simulation cell
const void Media::generate_simple(const Medium::Sizes& ar_sizes, const double z) {
    // Reserve memory for atoms
    reserve(4);

    // Add 4 atoms to the corners of simulation cell
    // The insertion order determines the orientation of big surface triangles
    add_atom( Point3(ar_sizes.xmin, ar_sizes.ymin, z) );
    add_atom( Point3(ar_sizes.xmax, ar_sizes.ymin, z) );
    add_atom( Point3(ar_sizes.xmax, ar_sizes.ymax, z) );
    add_atom( Point3(ar_sizes.xmin, ar_sizes.ymax, z) );

    calc_statistics();
}

// Generate edge with regular atom distribution between surface corners
const void Media::generate_middle(const Medium::Sizes& ar_sizes, const double z, const double r_cut) {
    const int n_atoms_per_side_x = ar_sizes.xbox / r_cut + 1;
    const int n_atoms_per_side_y = ar_sizes.ybox / r_cut + 1;

    // Reserve memory for atoms
    reserve( 2 * (n_atoms_per_side_x + n_atoms_per_side_y) - 4 );

    // Add atoms in x-edge
    for (int i = 1; i < n_atoms_per_side_x; ++i) {
        double y = ar_sizes.xmin + i * ar_sizes.xbox / n_atoms_per_side_x;
        add_atom( Point3(ar_sizes.xmin, y, z) );
        add_atom( Point3(ar_sizes.xmax, y, z) );
    }

    // Add atoms in y-edge
    for (int i = 1; i < n_atoms_per_side_y; ++i) {
        double x = ar_sizes.ymin + i * ar_sizes.ybox / n_atoms_per_side_y;
        add_atom( Point3(x, ar_sizes.ymin, z) );
        add_atom( Point3(x, ar_sizes.ymax, z) );
    }
}

// Extract surface by the atom types
const void Media::extract(const AtomReader& reader, const int type, const bool invert) {
    const int n_atoms = reader.get_n_atoms();
    vector<bool> is_type(n_atoms);

    // Get number and locations of atoms of desired type
    if (!invert)
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = reader.get_type(i) == type;
    else
        for (int i = 0; i < n_atoms; ++i)
            is_type[i] = reader.get_type(i) != type;

    // Preallocate memory for atoms
    reserve(vector_sum(is_type));

    // Store the atoms
    for (int i = 0; i < n_atoms; ++i)
        if (is_type[i])
            add_atom(reader.get_atom(i));
            
    calc_statistics();        
}

// Function extend the flat area by generating additional atoms
const Media Media::stretch(const double radius, const double box_width) {
    const double PI = 3.141592653589793;
    const int n_atoms = get_n_atoms();
    const double radius2 = radius * radius;
    const double box_w = (box_width/2)*sizes.zbox;
    const double radius_in = min(sizes.xbox/2.0, sizes.ybox/2.0);
        
    const int n_theta = 6;
    const double dtheta = 2*PI / n_theta;
    const int n_radius = (int) box_w * box_w / (1 * n_theta);
    
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);

    Media stretched(n_atoms + n_radius * n_theta);
    for (int i = 0; i < n_atoms; ++i)
        stretched.add_atom(get_point(i));
    
    // if input surface already is sufficiently wide, don't modify system at all
    if (box_w <= radius_in)
        return stretched;
    
    Media flat(n_atoms);  
    for (int i = 0; i < n_atoms; ++i)
        if (origin2d.distance2(get_point2(i)) > radius2)
            flat.add_atom(get_point(i));    
    
    flat.calc_statistics();
    Vec3 r0(sizes.xmid, sizes.ymid, flat.sizes.zmean);

    const double xmin = r0.x - box_w;
    const double xmax = r0.x + box_w;
    const double ymin = r0.y - box_w;
    const double ymax = r0.y + box_w;
   
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> radius_generator(0, 1);
    
    for (int j = 0; j < n_theta; ++j) {
        std::uniform_real_distribution<> theta_generator(j*dtheta, (j+1)*dtheta);
        for (int i = 0; i < n_radius; ++i) {      
            double r = radius_generator(gen) * (sqrt(2) * box_w - radius_in) + radius_in;        
            double theta = theta_generator(gen);
            
            double x = r * cos(theta) + origin2d.x;
            double y = r * sin(theta) + origin2d.y;
            
            bool in_box = x >= xmin && x <= xmax && y >= ymin && y <= ymax;
            bool in_flat = x >= sizes.xmin && x <= sizes.xmax && y >= sizes.ymin && y <= sizes.ymax;
            
            if (in_box && !in_flat)
                stretched.add_atom( Point3(x, y, flat.sizes.zmean) );
        }
    }
 
    stretched.calc_statistics();

    return stretched;
}

// Function to stretch the flat area
const Media Media::stretch_by_stretch(const double radius, const double coarse_factor) {
    const int n_atoms = get_n_atoms();
    const double radius2 = radius * radius;
    const double box_xy = 200.0;
    
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);

    Media stretched(n_atoms);
    Media flat(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        if (origin2d.distance2(get_point2(i)) <= radius2)
            stretched.add_atom(get_point(i));
        else
            flat.add_atom(get_point(i));
    }

    flat.calc_statistics();
    Vec3 r0(sizes.xmid, sizes.ymid, flat.sizes.zmean);

    const double xmin = r0.x - box_xy;
    const double xmax = r0.x + box_xy;
    const double ymin = r0.y - box_xy;
    const double ymax = r0.y + box_xy;

    for (int i = 0; i < flat.get_n_atoms(); ++i) {
        Point3 point = flat.get_point(i);
        Vec3 r(point.x, point.y, point.z);
        Vec3 dr = r - r0;
        double distance = dr.norm() - radius;
        double w = coarse_factor*pow(distance, 1.0);

//        dr = dr.normalize();
        dr *= w;
        point += Point3(dr.x, dr.y, dr.z);

//        double sign_x = 1.0;
//        if (dr.x < 0) sign_x = -1.0;
//        double sign_y = 1.0;
//        if (dr.y < 0) sign_y = -1.0;
//        point += Point3(w*sign_x, w*sign_y, 0);

        if (point.x >= xmin && point.x <= xmax && point.y >= ymin && point.y <= ymax)
            stretched.add_atom(point);
    }

    stretched.calc_statistics();

    return stretched;
}

// Function to coarsen the atoms with coarsener
const Media Media::coarsen(Coarseners &coarseners, const Medium::Sizes& ar_sizes) {
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);

    Media corners; corners.generate_simple(ar_sizes, coarseners.zmean);
    Media middle; middle.generate_middle(ar_sizes, coarseners.zmean, coarseners.r0_inf);
    middle.sort_atoms(3, "down", origin2d);
    this->sort_atoms(3, 2, "down", origin2d);

    Media union_surf;
    union_surf += corners;
    union_surf += middle;
    union_surf.add(this);

    return union_surf.clean(coarseners);
}

// Function to flatten the atoms on the sides of simulation box
const Media Media::rectangularize(const AtomReader::Sizes& ar_sizes, const double eps, const double latconst) {
    Coarseners coarseners;
    coarseners.attach_coarsener(make_shared<ConstCoarsener>(latconst / 3.0));

    Edge edge;
    edge.extract(this, ar_sizes, eps);

    Media surf;
    surf += edge;
    surf.add(this);

    return surf.clean(coarseners);
}

// Clean the surface from atoms that are too close to each other
const Media Media::clean(Coarseners &coarseners) {
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

    Media surf( n_atoms - vector_sum(do_delete) );
    for (int i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            surf.add_atom(get_atom(i));

    surf.calc_statistics();

    return surf;
}

// Function to delete atoms that are separate from others
// Atom is considered lonely if its coordination is lower than coord_min
const Media Media::clean_lonely_atoms(const double r_cut) {
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

    Media surf( vector_sum(atom_not_lonely) );
    for (int i = 0; i < n_atoms; ++i)
        if (atom_not_lonely[i])
            surf.add_atom(get_atom(i));

    surf.calc_statistics();
    return surf;
}

inline double Media::smooth_function(double distance, double smooth_factor) const {
    const double a = 1.0;
    return a * exp(-1.0 * distance / smooth_factor);
}

const void Media::smoothen(double radius, double smooth_factor, double r_cut) {
    // Calculate the horizontal span of the surface
    calc_statistics();
    Point2 origin2d(sizes.xmid, sizes.ymid);
    smoothen(origin2d, radius, smooth_factor, r_cut);
}

const void Media::smoothen(const Point2 &origin, double radius, double smooth_factor, double r_cut) {
    if (smooth_factor < 0.01) return;

    const int n_atoms = get_n_atoms();
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
            temp_surf.add_atom(get_point(i));

    // Smoothen the surface
    temp_surf.smoothen(smooth_factor, r_cut);

    // Write smoothened points back to surface
    int j = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (hot_atom[i])
            set_point(i, temp_surf.get_point(j++));
}

const void Media::smoothen(double smooth_factor, double r_cut) {
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
const void Edge::extract(const Medium* atoms, const AtomReader::Sizes& ar_sizes, const double eps) {
    const int n_atoms = atoms->get_n_atoms();

    // Reserve memory for atoms
    reserve(n_atoms);

    // Get the atoms from edge areas
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms->get_point(i);
        const bool near1 = on_boundary(point.x, ar_sizes.xmin, ar_sizes.xmax, eps);
        const bool near2 = on_boundary(point.y, ar_sizes.ymin, ar_sizes.ymax, eps);

        if (near1 || near2)
            add_atom(atoms->get_atom(i));
    }

    // Loop through all the added atoms and flatten the atoms on the sides of simulation box
    for (int i = 0; i < get_n_atoms(); ++i) {
        Point3 point = get_point(i);
        if ( on_boundary(point.x, ar_sizes.xmin, eps) ) set_x(i, ar_sizes.xmin);
        if ( on_boundary(point.x, ar_sizes.xmax, eps) ) set_x(i, ar_sizes.xmax);
        if ( on_boundary(point.y, ar_sizes.ymin, eps) ) set_y(i, ar_sizes.ymin);
        if ( on_boundary(point.y, ar_sizes.ymax, eps) ) set_y(i, ar_sizes.ymax);
    }
}

} /* namespace femocs */
