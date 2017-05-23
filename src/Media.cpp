/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Media.h"
#include "VoronoiMesh.h"
#include <numeric>

using namespace std;
namespace femocs {

Media::Media() : Medium() {}

Media::Media(const int n_atoms) : Medium(n_atoms) {}

Media::Media(const Medium::Sizes& s, const double z) {
    generate_simple(s, z);
}

// Generate system with 4 atoms at the corners
void Media::generate_simple(const Medium::Sizes& s, const double z) {
    // Reserve memory for atoms
    reserve(4);

    // Add 4 atoms to the corners of simulation cell
    // The insertion order determines the orientation of big surface triangles
    append( Point3(s.xmin, s.ymin, z) );
    append( Point3(s.xmax, s.ymin, z) );
    append( Point3(s.xmax, s.ymax, z) );
    append( Point3(s.xmin, s.ymax, z) );

    calc_statistics();
}

// Generate edge with regular atom distribution between surface corners
void Media::generate_middle(const Medium::Sizes& s, const double z, const double dist) {
    require(dist > 0, "Invalid distance between atoms: " + to_string(dist));
    const int n_atoms_per_side_x = s.xbox / dist + 1;
    const int n_atoms_per_side_y = s.ybox / dist + 1;
    const double dx = s.xbox / n_atoms_per_side_x;
    const double dy = s.ybox / n_atoms_per_side_y;

    // Reserve memory for atoms
    reserve( 2 * (n_atoms_per_side_x + n_atoms_per_side_y) - 4 );

    // Add atoms along x-edge
    for (int i = 1; i < n_atoms_per_side_x; ++i) {
        double x = s.xmin + i * dx;
        append( Point3(x, s.ymin, z) );
        append( Point3(x, s.ymax, z) );
    }

    // Add atoms along y-edge
    for (int i = 1; i < n_atoms_per_side_y; ++i) {
        double y = s.ymin + i * dy;
        append( Point3(s.xmin, y, z) );
        append( Point3(s.xmax, y, z) );
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

// Extend the flat area by reading additional atoms
Media Media::extend(const string &file_name, Coarseners &coarseners) {
    AtomReader reader;
    reader.import_file(file_name);

    Media stretched(reader.size());
    stretched += reader;
    stretched.calc_statistics();
    stretched.sizes.zmean = stretched.sizes.zmin;

    return stretched.coarsen(coarseners);
}

// Extend the flat area by generating additional atoms
Media Media::extend(const double latconst, const double box_width, Coarseners &coarseners) {
    calc_statistics();
    const double desired_box_width = box_width * sizes.zbox;

    // Over estimate the number of generated points and reserve memory for them
    int n_generated = pow(desired_box_width / latconst + 1, 2);
    n_generated -= (sizes.xbox / latconst - 1) * (sizes.ybox / latconst - 1);
    Media stretched( max(0, n_generated) );

    const double Wx = (desired_box_width - sizes.xbox) / 2.0;  // generation area width
    const double Wy = (desired_box_width - sizes.ybox) / 2.0;

    // if the input surface isn't sufficiently wide, add atoms to it, otherwise just add the boundary nodes
    if (Wx > 0 && Wy > 0) {
        // add points outside the already existing area
        for (double y = sizes.ymin - Wy; y <= sizes.ymax + Wy; y += latconst)
            for (double x = sizes.xmin - Wx; x <= sizes.xmax + Wx; x += latconst)
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

// Coarsen the atoms by generating additional boundary nodes and then running cleaner
Media Media::coarsen(Coarseners &coarseners) {
    Media corners, middle, union_surf;

    corners.generate_simple(sizes, sizes.zmean);
    middle.generate_middle( sizes, sizes.zmean, coarseners.get_r0_inf(sizes) );
    sort_atoms(3, "down");

    union_surf += corners;
    union_surf += middle;
    union_surf += *this;

    return union_surf.clean(coarseners);
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
void Media::faces_clean(const TetgenMesh& mesh, const double r_cut) {
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

// Extract the surface atoms whose Voronoi cells are exposed to vacuum
bool Media::voronoi_clean(vector<Vec3>& areas, const double radius, const double latconst, const string& mesh_quality) {
    const int n_atoms = size();

    // Separate nanotip from the whole surface
    Media nanotip;
    get_nanotip(nanotip, radius);

    // Generate Voronoi cells around the nanotip
    VoronoiMesh voromesh;
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    if( voromesh.generate(nanotip, latconst, "rQq" + mesh_quality, "vQ") )
        return 1;
    
    // Clean the mesh from faces and cell that have node in the infinity
    voromesh.clean();
    
    // Extract the surface faces and cells
    voromesh.mark_mesh(nanotip);
    
    voromesh.nodes.write("output/voro_nodes.vtk");
    voromesh.vfaces.write("output/voro_faces.vtk");
    voromesh.voros.write("output/voro_cells.vtk");
    
    // Get the atoms and their surface areas whose Voronoi cells are exposed to vacuum
    voromesh.extract_surface(*this, areas, nanotip);
    calc_statistics();
    
    return 0;
}

// Separate cylindrical region from substrate region
void Media::get_nanotip(Media& nanotip, const double radius) {
    const int n_atoms = size();
    const double radius2 = radius * radius;

    // Make map for atoms in nanotip
    vector<bool> is_nanotip; is_nanotip.reserve(n_atoms);
    Point2 centre(sizes.xmid, sizes.ymid);

    for (int i = 0; i < n_atoms; ++i)
       is_nanotip.push_back( centre.distance2(get_point2(i)) <= radius2 );

    // Reserve memory for nanotip and substrate
    const int n_nanotip = vector_sum(is_nanotip);
    nanotip.reserve(n_nanotip);
    vector<Atom> atoms_save; atoms_save.reserve(n_atoms - n_nanotip);

    // Separate nanotip and substrate
    for (int i = 0; i < n_atoms; ++i) {
        if (is_nanotip[i])
            nanotip.append(get_atom(i));
        else
            atoms_save.push_back(get_atom(i));
    }

    nanotip.calc_statistics();
    atoms = atoms_save;
}

// Smoothen the atoms inside the cylinder
void Media::smoothen(const double radius, const double smooth_factor, const double r_cut) {
    if (smooth_factor <= 0) return;

    // Calculate the horizontal span of the surface
    calc_statistics();

    Media nanotip;
    get_nanotip(nanotip, radius);
    nanotip.smoothen(smooth_factor, r_cut);

    *this += nanotip;
}

// Smoothen all the atoms in the system
void Media::smoothen(const double smooth_factor, const double r_cut) {
    if (smooth_factor <= 0) return;

    const double r_cut2 = r_cut * r_cut;
    const double decay_factor = -1.0 / smooth_factor;
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

            double weight = exp(decay_factor * sqrt(distance2));
            atoms[i].point += point2 * weight;
            atoms[j].point += point1 * weight;
            weights_sum[i] += weight;
            weights_sum[j] += weight;
        }
    }

    // Normalise smoothed vertices
    for (int i = 0; i < n_atoms; ++i) {
        if (weights_sum[i] > 0)
            atoms[i].point *= 1.0 / weights_sum[i];
        else
            atoms[i].point = points[i];
    }
}

} /* namespace femocs */
