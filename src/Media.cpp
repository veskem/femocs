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

void Media::transform(const double latconst) {
    const int n_atoms = size();
    calc_statistics();
    Point3 origin(sizes.xmid, sizes.ymid, sizes.zmid);

    double fx = 1.0 + 3.0 * latconst / sizes.xbox;
    double fy = 1.0 + 3.0 * latconst / sizes.ybox;
    double fz = 1.0 + 3.0 * latconst / sizes.zbox;
    Point3 df(fx, fy, fz);

    for (int i = 0; i < n_atoms; ++i)
        atoms[i].point *= df;

    calc_statistics();
    origin -= Point3(sizes.xmid, sizes.ymid, sizes.zmid);

    for (int i = 0; i < n_atoms; ++i)
        atoms[i].point += origin;
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
void Media::clean_by_triangles(vector<int>& surf2face, const TriangleInterpolator& interpolator, const double r_cut) {
    if (r_cut <= 0) return;

    const int n_atoms = size();
    vector<Atom> _atoms;
    _atoms.reserve(n_atoms);
    surf2face.clear();
    surf2face.reserve(n_atoms);

    for (int i = 0; i < n_atoms; ++i) {
        Atom atom = get_atom(i);
        atom.marker = interpolator.near_surface(atom.point, r_cut);
        if (atom.marker >= 0) {
            _atoms.push_back(atom);
            surf2face.push_back(atom.marker);
        }
    }

    atoms = _atoms;
    calc_statistics();
}

int Media::calc_voronois(VoronoiMesh& voromesh, vector<bool>& node_in_nanotip,
        const double radius, const double latconst, const string& mesh_quality)
{
    const int n_this_nodes = size();
    const double radius2 = radius * radius;
    Medium::calc_statistics();

    // Make map for atoms in nanotip
    Point2 centre(sizes.xmid, sizes.ymid);
    node_in_nanotip = vector<bool>(n_this_nodes);
    for (int i = 0; i < n_this_nodes; ++i)
        node_in_nanotip[i] = centre.distance2(get_point2(i)) <= radius2;

    const int n_nanotip_nodes = vector_sum(node_in_nanotip);

    // Separate nanotip from substrate
    Medium nanotip(n_nanotip_nodes);
    for (int i = 0; i < n_this_nodes; ++i)
        if (node_in_nanotip[i])
            nanotip.append(get_atom(i));

    nanotip.calc_statistics();

    double t0;
    start_msg(t0, "  Generating Voronoi mesh...");

    // Generate Voronoi cells around the nanotip
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    int errcode = voromesh.generate(nanotip, latconst, "rQq" + mesh_quality, "vQ");
    if (errcode) return errcode;

    // Clean the mesh from faces and cells that have node in the infinity
    voromesh.clean();
    end_msg(t0);

    voromesh.nodes.write("out/voro_nodes.vtk");
    voromesh.vfaces.write("out/voro_faces.vtk");
    voromesh.voros.write("out/voro_cells.vtk");

    return 0;
}

int Media::clean_by_voronois(const double radius, const double latconst, const string& mesh_quality) {
    // Extract nanotip
    Media nanotip;
    vector<bool> node_in_nanotip;
    const int n_nanotip_nodes = get_nanotip(nanotip, node_in_nanotip, radius);

    // Generate Voronoi cells around the nanotip
    VoronoiMesh mesh;

    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    // F - suppress output of faces and edges, B - suppress output of boundary info
    const int err_code = mesh.generate(nanotip, latconst, "rQFBq" + mesh_quality, "vQFB");
    if (err_code) return err_code;

    // Clean the mesh from faces and cells that have node in the infinity
    mesh.clean();

    // Extract the surface faces and cells
    mesh.mark_mesh(nanotip, latconst);

//    mesh.nodes.write("out/voro_nodes.vtk");
//    mesh.vfaces.write("out/voro_faces.vtk");
//    mesh.voros.write("out/voro_cells.vtk");

    require(mesh.nodes.size() > 0, "Empty Voronoi mesh cannot be handled!");
    require(mesh.voros.size() > 0, "Empty Voronoi mesh cannot be handled!");

    // delete atoms whose Voronoi cell is not exposed to vacuum
    vector<Atom> _atoms;
    _atoms.reserve(size());

    int cell = 0;
    for (int i = 0; i < size(); ++i)
        if (node_in_nanotip[i]) {
            if (mesh.voros.get_marker(cell++) == TYPES.SURFACE)
                _atoms.push_back(get_atom(i));
        } else
            _atoms.push_back(get_atom(i));

    atoms = _atoms;
    calc_statistics();

    return 0;
}

// Separate cylindrical region from substrate region
int Media::get_nanotip(Media& nanotip, vector<bool>& node_in_nanotip, const double radius) {
    const int n_this_nodes = size();
    const double radius2 = radius * radius;
    Medium::calc_statistics();

    // Make map for atoms in nanotip
    Point2 centre(sizes.xmid, sizes.ymid);
    node_in_nanotip = vector<bool>(n_this_nodes);
    for (int i = 0; i < n_this_nodes; ++i)
        node_in_nanotip[i] = centre.distance2(get_point2(i)) <= radius2;

    const int n_nanotip_nodes = vector_sum(node_in_nanotip);

    // Separate nanotip from substrate
    nanotip.reserve(n_nanotip_nodes);
    for (int i = 0; i < n_this_nodes; ++i)
        if (node_in_nanotip[i])
            nanotip.append(get_atom(i));

    nanotip.calc_statistics();
    return n_nanotip_nodes;
}

// Extract the surface atoms whose Voronoi cells are exposed to vacuum
int Media::voronoi_clean(vector<Vec3>& areas, const double radius, const double latconst, const string& mesh_quality) {
    // Separate nanotip from the whole surface
    Media nanotip;
    get_nanotip(nanotip, radius);

    // Generate Voronoi cells around the nanotip
    VoronoiMesh voromesh;
    // r - reconstruct, v - output Voronoi cells, Q - quiet, q - mesh quality
    // F - suppress output of faces and edges, B - suppress output of boundary info
    const int err_code = voromesh.generate(nanotip, latconst, "rQFBq" + mesh_quality, "vQFB");
    if(err_code) return err_code;

    // Clean the mesh from faces and cell that have node in the infinity
    voromesh.clean();

    // Extract the surface faces and cells
    voromesh.mark_mesh(nanotip, latconst);

    voromesh.nodes.write("out/voro_nodes.vtk");
    voromesh.vfaces.write("out/voro_faces.vtk");
    voromesh.voros.write("out/voro_cells.vtk");

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

// Smoothen the atoms inside the cylinder
void Media::smoothen(const Config& conf, const double r_cut) {
    if (r_cut <= 0) return;

    // Calculate the horizontal span of the surface
    calc_statistics();

    Media nanotip;
    get_nanotip(nanotip, conf.radius);

    vector<vector<unsigned>> nborlist;
    nanotip.calc_nborlist(nborlist, conf.nnn, r_cut);

    for (int i = 0; i < 3; ++i) {
        nanotip.laplace_smooth(conf.smooth_lambda, nborlist);
        nanotip.laplace_smooth(conf.smooth_mu, nborlist);
    }

    *this += nanotip;
}

// Apply one cycle of Taubin lambda|mu algorithm
void Media::laplace_smooth(const double scale, const vector<vector<unsigned>>& nborlist) {
    size_t n_nodes = size();
    vector<Point3> displacements(n_nodes);

    // Get per-vertex displacement
    for (size_t i = 0; i < n_nodes; ++i) {
        // Skip lonely vertices
        if (nborlist[i].size() == 0)
            continue;

        const double weight = 1.0 / nborlist[i].size();

        // Sum the displacements
        Point3 point = get_point(i);
        for(size_t nbor : nborlist[i])
            displacements[i] += (get_point(nbor) - point) * weight;
    }

    // Apply per-point displacement
    for (size_t i = 0; i < n_nodes; ++i)
        atoms[i].point += displacements[i] * scale;
}


// Calculate list of close neighbours using brute force technique
void Media::calc_nborlist(vector<vector<unsigned>>& nborlist, const int nnn, const double r_cut) {
    require(r_cut > 0, "Invalid cut-off radius: " + to_string(r_cut));

    const size_t n_atoms = size();
    const double r_cut2 = r_cut * r_cut;
    const double eps = 0.001 * r_cut;

    // Initialise list of closest neighbours
    nborlist = vector<vector<unsigned>>(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        nborlist[i].reserve(nnn);

    // Loop through all the atoms
    for (size_t i = 0; i < n_atoms - 1; ++i) {

        Point3 point1 = get_point(i);

        // Skip the points that are on the boundary of simubox
        if (on_boundary(point1.x, sizes.xmin, sizes.xmax, eps) ||
                on_boundary(point1.y, sizes.ymin, sizes.ymax, eps))
            continue;

        // Loop through all the possible neighbours of the atom
        for (size_t j = i + 1; j < n_atoms; ++j) {
            if ( r_cut2 >= point1.distance2(get_point(j)) ) {
                nborlist[i].push_back(j);
                nborlist[j].push_back(i);
            }
        }
    }
}
} /* namespace femocs */
