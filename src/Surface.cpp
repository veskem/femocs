/*
 * Media.cpp
 *
 *  Created on: 5.4.2016
 *      Author: veske
 */

#include "Surface.h"
#include "AtomReader.h"
#include <numeric>

using namespace std;
namespace femocs {

Surface::Surface() : Medium(), coarseners(NULL) {}

Surface::Surface(const int n_atoms) : Medium(n_atoms), coarseners(NULL) {}

Surface::Surface(const Medium::Sizes& s, const double z) :
        Medium(4), coarseners(NULL)
{
    // Add 4 atoms to the corners of simulation cell
    // The insertion order determines the orientation of big surface triangles
    append( Point3(s.xmin, s.ymin, z) );
    append( Point3(s.xmax, s.ymin, z) );
    append( Point3(s.xmax, s.ymax, z) );
    append( Point3(s.xmin, s.ymax, z) );

    calc_statistics();
}

Surface::Surface(const Medium::Sizes& s, const double z, const double dist) :
        coarseners(NULL)
{
    require(dist > 0, "Invalid distance between atoms: " + d2s(dist));

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

void Surface::extend(Surface& extension, double latconst,
    double box_width, double z, const Sizes& sizes, bool circular)
{
    require(coarseners, "Coarseners not provided!");
    require(coarseners->get_r0_inf(sizes) > 0, "Invalid coarsener detected.");

    const double desired_box_width = box_width * sizes.zbox;

    // if the input surface isn't sufficiently wide, add atoms to it
    if (desired_box_width > sizes.xbox && desired_box_width > sizes.ybox) {
        // Over estimate the number of generated points and reserve memory for them
        int n_generated = pow(desired_box_width / latconst + 1, 2);
        n_generated -= (sizes.xbox / latconst - 1) * (sizes.ybox / latconst - 1);
        extension.reserve(max(0, n_generated));

        const double r_max = sqrt(2) * desired_box_width / 2.0;
        double r_min;
        if (circular) r_min = coarseners->get_radius();
        else r_min = min(sizes.xbox, sizes.ybox) / 2.0;

        double dR = coarseners->get_cutoff(Point3(sizes.xmid - r_min, sizes.ymid - r_min, z));
        int i;

        // add points along a rectangular spiral
        for (double r = r_min; r <= r_max; r += dR) {
            const int n_points = int(2.0 * r / dR);
            dR = 2.0 * r / n_points;

            double x = sizes.xmid - r;
            double y = sizes.ymid - r;

            for (i = 0; i < n_points; ++i, y += dR)
                extension.append(Point3(x, y, z));
            for (i = 0; i < n_points; ++i, x += dR)
                extension.append(Point3(x, y, z));
            for (i = 0; i < n_points; ++i, y -= dR)
                extension.append(Point3(x, y, z));
            for (i = 0; i < n_points; ++i, x -= dR)
                extension.append(Point3(x, y, z));

            dR = coarseners->get_cutoff(Point3(x - dR, y - dR, z));
        }
    }

    // if the input already is sufficiently wide, add just the boundary nodes
    else {
        Surface middle(sizes, z, coarseners->get_r0_inf(sizes));
        extension = Surface(sizes, z);
        extension += middle;
    }

    extension.calc_statistics();
}

void Surface::coarsen(Surface& surface) {
    require(coarseners, "Coarseners not provided!");
    sort_atoms(3, "down");

    Surface middle(sizes, sizes.zmean, coarseners->get_r0_inf(sizes));
    Surface union_surf(sizes, sizes.zmean);
    union_surf += middle;
    union_surf += surface;

    clean(union_surf);
}

void Surface::clean(Surface& surface) {
    require(coarseners, "Coarseners not provided!");
    const int n_atoms = surface.size();
    vector<bool> do_delete(n_atoms, false);

    // Loop through all the atoms
    for (int i = 0; i < n_atoms - 1; ++i) {
        // Skip already deleted atoms
        if (do_delete[i]) continue;

        Point3 point1 = surface.get_point(i);
        coarseners->pick_cutoff(point1);

        for (int j = i+1; j < n_atoms; ++j) {
            // Skip already deleted atoms
            if (do_delete[j]) continue;
            do_delete[j] = coarseners->nearby(point1, surface.get_point(j));
        }
    }

    // remove coarsened atoms
    int j = 0;
    for (int i = 0; i < n_atoms; ++i)
        if(!do_delete[i])
            surface.atoms[j++] = surface.get_atom(i);

    surface.atoms.resize(j);
    surface.calc_statistics();
}

void Surface::add_cleaned_roi_to(Surface& surface) {
    require(coarseners, "Coarseners not provided!");
    const int n_atoms = size();
    vector<int> do_delete(n_atoms, 0);

    // mark atoms outside the nanotip
    for (int i = 0; i < n_atoms; ++i)
        do_delete[i] = -1 * !coarseners->inside_roi(get_point(i));

    // Loop through all the nanotip atoms
    for (int i = 0; i < n_atoms; ++i) {
        // skip already marked atoms
        if (do_delete[i] != 0) continue;

        Point3 point1 = get_point(i);
        coarseners->pick_cutoff(point1);

        for (int j = i+1; j < n_atoms; ++j) {
            // skip already marked atoms
            if (do_delete[j] == 0)
                do_delete[j] = coarseners->nearby(point1, get_point(j));
        }
    }

    // add coarsened atoms to the input surface

    int n_coarsened_atoms = 0;
    for (int dd : do_delete)
        if (dd <= 0)
            n_coarsened_atoms++;

    surface.atoms.reserve(surface.size() + n_coarsened_atoms);
    for (int i = 0; i < n_atoms; ++i)
        if(do_delete[i] <= 0)
            surface.append(get_atom(i));

    surface.calc_statistics();
}

/* TODO: Leaves bigger holes into system than brute force method,
 * because the atoms in linked list are not radially ordered. Do something about it! */
void Surface::fast_coarsen(Surface &surface, const Medium::Sizes &s) {
    require(coarseners, "Coarseners not provided!");
    calc_linked_list(coarseners->get_r0_inf(s));

    const int n_atoms = size();
    require((int)list.size() == n_atoms, "Invalid linked list size: " + d2s(list.size()));
    require((int)head.size() == nborbox_size[0]*nborbox_size[1]*nborbox_size[2],
            "Invalid linked list header size: " + d2s(head.size()));

    vector<int> do_delete(n_atoms, false);

    // loop through the atoms
    for (int i = 0; i < n_atoms; ++i) {
        // skip the atoms that are already deleted
        if (do_delete[i]) continue;

        array<int,3>& i_atom = nborbox_indices[i];

        Point3 point1 = atoms[i].point;
        coarseners->pick_cutoff(point1);

        // loop through the boxes where the neighbours are located; there are up to 3^3=27 boxes
        for (int iz = i_atom[2]-1; iz <= i_atom[2]+1; ++iz) {
            // some of the iterations are be skipped if the box is on the simu box boundary
            if (iz < 0 || iz >= nborbox_size[2]) continue;
            for (int iy = i_atom[1]-1; iy <= i_atom[1]+1; ++iy) {
                if (iy < 0 || iy >= nborbox_size[1]) continue;
                for (int ix = i_atom[0]-1; ix <= i_atom[0]+1; ++ix) {
                    if (ix < 0 || ix >= nborbox_size[0]) continue;

                    // transform volumetric neighbour box index to linear one
                    int i_cell = (iz * nborbox_size[1] + iy) * nborbox_size[0] + ix;
                    require(i_cell >= 0 && i_cell < (int)head.size(),
                            "Invalid neighbouring cell index: " + d2s(i_cell));

                    // get the index of first atom in given neighbouring cell
                    int j = head[i_cell];

                    // loop through all atoms in a given neighbouring cell
                    while(j >= 0) {
                        require(j < n_atoms, "Invalid index in linked list: " + d2s(j));
                        // skip the same atoms and the atoms that are already deleted
                        if (!do_delete[j] && i != j)
                            do_delete[j] = coarseners->nearby(point1, get_point(j));
                        j = list[j];
                    }
                }
            }
        }
    }

    // Store coarsened surface
    surface.reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        if (!do_delete[i])
            surface.append(get_atom(i));
    surface.calc_statistics();
}

void Surface::clean_by_triangles(Interpolator& interpolator, const TetgenMesh* mesh, const double r_cut) {
    if (r_cut <= 0) return;

    const int n_atoms = size();

    interpolator.lintri.set_mesh(mesh);
    interpolator.lintri.precompute();

    int face = 0, j = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Atom atom = get_atom(i);
        face = abs(interpolator.lintri.locate_cell(atom.point, face));
        if (interpolator.lintri.fast_distance(atom.point, face) < r_cut) {
            atom.marker = face;
            atoms[j++] = atom;
        }
    }

    atoms.resize(j);
    calc_statistics();
}

void Surface::get_nanotip(Surface& nanotip, const Coarseners& coarseners) {
    const int n_atoms = size();

    // Make map for atoms in nanotip
    vector<bool> is_nanotip(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
       is_nanotip[i] = coarseners.inside_roi(get_point(i));

    // Reserve memory for nanotip and substrate
    nanotip.reserve(vector_sum(is_nanotip));

    // Separate nanotip and substrate
    int j = 0;
    for (int i = 0; i < n_atoms; ++i) {
        if (is_nanotip[i])
            nanotip.append(get_atom(i));
        else
            atoms[j++] = get_atom(i);
    }

    atoms.resize(j);
    nanotip.calc_statistics();
}

void Surface::smoothen_roi(double smooth_factor, double r_cut) {
    if (smooth_factor <= 0) return;
    require(coarseners, "Coarseners not provided!");

    // Calculate the horizontal span of the surface
    calc_statistics();

    Surface nanotip;
    get_nanotip(nanotip, *coarseners);
    nanotip.smoothen(smooth_factor, r_cut);

    *this += nanotip;
}

void Surface::smoothen(double smooth_factor, double r_cut) {
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

void Surface::extend(Surface& extended_surf, const Config& conf) {
    require(coarseners, "Coarseners not provided!");

    if (conf.path.extended_atoms == "") {
        // Extend surface by generating additional nodes
        calc_statistics();
        extend(extended_surf, conf.geometry.latconst, conf.geometry.box_width,
            coarseners->centre.z, sizes, false);
    }

    else {
        // Extend surface by first reading the data from file...
        AtomReader reader(&conf.geometry);
        reader.import_file(conf.path.extended_atoms);
        extended_surf += reader;
        extended_surf.sort_atoms(3, "up");

        // ... and then by adding points on horizontal plane, if necessary
        Surface temp_surf;
        extend(temp_surf, conf.geometry.latconst, conf.geometry.box_width,
            extended_surf.sizes.zmin, reader.sizes, true);
        extended_surf += temp_surf;
        clean(extended_surf);
        extended_surf.sizes.zmean = extended_surf.sizes.zmin;
    }
}

int Surface::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum,
        const Surface& extended_surf, const Config& conf)
{
    // sort atoms radially to increase the symmetry of the resulting surface
    sort_atoms(3, "up");

    // Coarsen & smoothen surface
    coarse_surf.set_coarsener(coarseners);
    coarse_surf.atoms = extended_surf.atoms;
//    coarse_surf += *this;
    add_cleaned_roi_to(coarse_surf);
    clean(coarse_surf);
    coarse_surf.smoothen_roi(conf.smoothing.beta_atoms, 3.0*conf.geometry.coordination_cutoff);

    // Generate bulk & vacuum corners
    coarse_surf.calc_statistics();  // calculate zmin and zmax for surface
    double vacuum_height = max(conf.geometry.latconst, coarse_surf.sizes.zbox) * conf.geometry.box_height;
    double bulk_height = conf.geometry.bulk_height * conf.geometry.latconst;

    vacuum = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin + vacuum_height);
    bulk   = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin - bulk_height);

    return 0;
}

} /* namespace femocs */
