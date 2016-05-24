/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

using namespace std;
namespace femocs {

Mesher::Mesher(const string mesher) {
    this->mesher = mesher;
}

// Function to generate simple mesh that consists of one tetrahedron
const void Mesher::get_test_mesh(Mesh* new_mesh) {
    const int n_nodes = 4;
    const int n_faces = 4;
    const int n_elems = 1;

    new_mesh->init_nodes(n_nodes);
    new_mesh->add_node(Point3d(1.0, 0.0, 0.7));
    new_mesh->add_node(Point3d(-1.0, 0.0, 0.7));
    new_mesh->add_node(Point3d(0.0, 1.0, -0.7));
    new_mesh->add_node(Point3d(0.0, -1.0, -0.7));

    new_mesh->init_faces(n_faces);
    new_mesh->add_face(0, 1, 3);
    new_mesh->add_face(1, 2, 3);
    new_mesh->add_face(2, 0, 3);
    new_mesh->add_face(0, 1, 2);

    new_mesh->init_elems(n_elems);
    new_mesh->add_elem(0, 1, 2, 3);
}

// Function to generate mesh from surface, bulk and vacuum atoms
const void Mesher::get_volume_mesh(Mesh* new_mesh, Bulk* bulk, Surface* surf, Vacuum* vacuum,
        const string cmd) {
    int i;
    int n_bulk = bulk->get_n_atoms();
    int n_surf = surf->get_n_atoms();
    int n_vacuum = vacuum->get_n_atoms();

    new_mesh->init_nodes(n_bulk + n_surf + n_vacuum);

    for (i = 0; i < n_surf; ++i)
        new_mesh->add_node(surf->get_point(i));

    for (i = 0; i < n_bulk; ++i)
        new_mesh->add_node(bulk->get_point(i));

    for (i = 0; i < n_vacuum; ++i)
        new_mesh->add_node(vacuum->get_point(i));

    new_mesh->indxs.surf_start = 0;
    new_mesh->indxs.surf_end = n_surf - 1;
    new_mesh->indxs.bulk_start = n_surf;
    new_mesh->indxs.bulk_end = n_surf + n_bulk - 1;
    new_mesh->indxs.vacuum_start = n_surf + n_bulk;
    new_mesh->indxs.vacuum_end = n_surf + n_bulk + n_vacuum - 1;
    new_mesh->indxs.tetgen_start = n_surf + n_bulk + n_vacuum;

    new_mesh->recalc(cmd);
}

const vector<bool> Mesher::get_vacuum_indices(Mesh* big_mesh, const int n_bulk, const int n_surf,
        const double zmin) {

    int i, j;
    int n_elems = big_mesh->get_n_elems();

    vector<bool> is_vacuum;
    is_vacuum.resize(n_elems);

    int n_vacuum, n_not_bottom, n_not_bulk, n_surface;

    for (i = 0; i < n_elems; ++i) {
        n_vacuum = n_not_bottom = n_not_bulk = n_surface = 0;

        // Find nodes that are not inside the bulk
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (big_mesh->get_elem(i, j) >= (n_bulk + n_surf)) n_vacuum++;

        for (j = 0; j < n_nodes_per_elem; ++j)
            if (big_mesh->get_elem(i, j) < n_surf) n_surface++;

        // Find nodes that are not in floating element
        for (j = 0; j < n_nodes_per_elem; ++j)
            if ((big_mesh->get_elem(i, j) >= (n_bulk + n_surf)) || (big_mesh->get_elem(i, j) < n_surf))
                n_not_bulk++;

//        // Find elements not belonging to the bottom of simulation cell
//        for (j = 0; j < n_nodes_per_elem; ++j) {
//            int node = big_mesh->get_elem(i, j);
//            if (fabs(big_mesh->get_z(node) - zmin) > 1.0 * latconst) n_not_bottom++;
//        }

        //is_vacuum[i] = (n_vacuum > 0) && (n_not_bulk > 1) && (n_not_bottom == n_nodes_per_elem);
        is_vacuum[i] = ((n_vacuum > 0)  || (n_surface > 3)) && (n_not_bulk > 1);// && (n_not_bulk > 1);
    }
    return is_vacuum;
}

const void Mesher::separate_meshes(Mesh* bulk, Mesh* vacuum, Mesh* big_mesh, const int n_bulk,
        const int n_surf, const double zmin, const string cmd) {

    int i, j;
    int n_nodes = big_mesh->get_n_nodes();
    int n_elems = big_mesh->get_n_elems();
    const int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum = get_vacuum_indices(big_mesh, n_bulk, n_surf, zmin);
    int n_vacuum_elems = vector_sum(elem_in_vacuum);

    // Copy the nodes from input mesh without modification
    vacuum->init_nodes(n_nodes);
    vacuum->copy_nodes(big_mesh);

    bulk->init_nodes(n_nodes);
    bulk->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_elems - n_vacuum_elems);

    for (i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i]) {
            j = n_nodes_per_elem * i;
            vacuum->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        } else {
            j = n_nodes_per_elem * i;
            bulk->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

const void Mesher::separate_meshes_bymarker(Mesh* bulk, Mesh* vacuum, Mesh* big_mesh,
        const AtomReader::Types* types, const string cmd) {

    int i, j;
    int n_nodes = big_mesh->get_n_nodes();
    int n_elems = big_mesh->get_n_elems();
    const int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum(n_elems);
    for (i = 0; i < n_elems; ++i)
        elem_in_vacuum[i] = big_mesh->get_elemmarker(i) != types->type_bulk;

    int n_vacuum_elems = vector_sum(elem_in_vacuum);

    // Copy the nodes from input mesh without modification
    vacuum->init_nodes(n_nodes);
    vacuum->copy_nodes(big_mesh);

    bulk->init_nodes(n_nodes);
    bulk->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_elems - n_vacuum_elems);

    for (i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i]) {
            j = n_nodes_per_elem * i;
            vacuum->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        } else {
            j = n_nodes_per_elem * i;
            bulk->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

// Function to manually generate surface faces into available mesh
const void Mesher::generate_monolayer_surf_faces(Mesh* mesh) {
    int i, j, n1, n2, n3;
    int locs[n_nodes_per_elem];
    bool b_locs[n_nodes_per_elem];

    int n_elems = mesh->get_n_elems();
    int n_nodes = mesh->get_n_nodes();
    int n_faces = mesh->get_n_faces();

    const int* elems = mesh->get_elems();
    vector<bool> is_surface(n_elems);

    // Mark the elements that have exactly one face on the surface and one node in the vacuum
    for (i = 0; i < n_elems; ++i) {
        for (j = 0; j < n_nodes_per_elem; ++j) {
            int elem = mesh->get_elem(i, j);
            // Node in bulk?
            if ((elem >= mesh->indxs.bulk_start) && (elem <= (mesh->indxs.bulk_end)))
                locs[j] = -1;

            // Node in vacuum?
            else if (elem >= mesh->indxs.vacuum_start)
                locs[j] = 0;

            // Node in surface!
            else
                locs[j] = 1;
        }

        // Bulk, surf and vacuum are encoded in such a way that
        // sum(locs) == 3 only if 3 nodes are on surface and 1 in vacuum
        is_surface[i] = locs[0] + locs[1] + locs[2] + locs[3] == 3;
    }

    int n_surf_faces = vector_sum(is_surface);

    // Make temporary copy of available faces
//    Mesh temp_mesh(this->mesher);
//    temp_mesh.init_faces(n_faces);
//    temp_mesh.copy_faces(mesh, 0);
//
//    // Make face list size bigger
//    mesh->init_faces(n_faces + n_surf_faces);
//    // Copy back already available faces
//    mesh->copy_faces(&temp_mesh, 0);

    mesh->init_faces(n_surf_faces);

    // Generate the faces that separate material and vacuum
    // It is assumed that each element has only one face on the surface
    for (i = 0; i < n_elems; ++i)
        if (is_surface[i]) {
            for (j = 0; j < n_nodes_per_elem; ++j)
                b_locs[j] = mesh->get_elem(i, j) <= mesh->indxs.bulk_end;

            /* The possible combinations of b_locs and n1,n2,n3:
             * b_locs: 1110  1101  1011  0111
             *     n1: elem0 elem0 elem0 elem1
             *     n2: elem1 elem1 elem2 elem2
             *     n3: elem2 elem3 elem3 elem3
             */
            j = n_nodes_per_elem * i;
            n1 = b_locs[0] * elems[j + 0] + (!b_locs[0]) * elems[j + 1];
            n2 = (b_locs[0] & b_locs[1]) * elems[j + 1] + (b_locs[2] & b_locs[3]) * elems[j + 2];
            n3 = (!b_locs[3]) * elems[j + 2] + b_locs[3] * elems[j + 3];
            mesh->add_face(n1, n2, n3);
        }
}

// Function to manually generate surface faces into available mesh
const void Mesher::generate_surf_faces(Mesh* mesh) {
    int elem, node, n0, n1, n2;

    const int n_elems = mesh->get_n_elems();
    const int max_surf_indx = mesh->indxs.surf_end;
    const int* elems = mesh->get_elems();
    vector<bool> elem_in_surface(n_elems);
    bool surf_locs[n_nodes_per_elem];

    // Mark the elements that have exactly one face on the surface
    for (elem = 0; elem < n_elems; ++elem) {
        int n_surf_nodes = 0;
        for (node = 0; node < n_nodes_per_elem; ++node)
            if (mesh->get_elem(elem, node) <= max_surf_indx)
                n_surf_nodes++;

        elem_in_surface[elem] = (n_surf_nodes == 3);
    }

    // Reserve memory for faces in the mesh
    mesh->init_faces(vector_sum(elem_in_surface));

    // Generate the faces that separate material and vacuum
    // It is assumed that each element has only one face on the surface
    for (elem = 0; elem < n_elems; ++elem)
        if (elem_in_surface[elem]) {
            // Find the indices of nodes that are on the surface
            for (node = 0; node < n_nodes_per_elem; ++node)
                surf_locs[node] = mesh->get_elem(elem, node) <= max_surf_indx;

            /* The possible combinations of surf_locs and n0,n2,n3:
             * surf_locs: 1110  1101  1011  0111
             *        n0: elem0 elem0 elem0 elem1
             *        n1: elem1 elem1 elem2 elem2
             *        n2: elem2 elem3 elem3 elem3
             */
            node = n_nodes_per_elem * elem;
            n0 = surf_locs[0] * elems[node + 0] + (!surf_locs[0]) * elems[node + 1];
            n1 = (surf_locs[0] & surf_locs[1]) * elems[node + 1] + (surf_locs[2] & surf_locs[3]) * elems[node + 2];
            n2 = (!surf_locs[3]) * elems[node + 2] + surf_locs[3] * elems[node + 3];
            mesh->add_face(n0, n1, n2);
        }
}

// Precompute the data needed to execute the Moller-Trumbore algorithm
const void Mesher::precompute_triangles(Mesh* mesh, const Vec3d &direction) {
    Vec3d v0, v1, v2, e1, e2, pv;
    int face, coord;
    double det, i_det;
    int n_faces = mesh->get_n_faces();

    precomp.edge1.reserve(n_faces);
    precomp.edge2.reserve(n_faces);
    precomp.vert0.reserve(n_faces);
    precomp.pvec.reserve(n_faces);
    precomp.is_parallel.reserve(n_faces);

    // Loop through all the faces
    for (face = 0; face < n_faces; ++face) {

        // Loop through all the coordinates
        for (coord = 0; coord < n_coordinates; ++coord) {
            v0[coord] = mesh->get_node(mesh->get_face(face, 0), coord);
            v1[coord] = mesh->get_node(mesh->get_face(face, 1), coord);
            v2[coord] = mesh->get_node(mesh->get_face(face, 2), coord);
        }

        e1 = v1 - v0;
        e2 = v2 - v0;
        pv = direction.crossProduct(e2);
        det = e1.dotProduct(pv);
        i_det = 1.0 / det;

        precomp.edge1[face] = e1 * i_det;
        precomp.edge2[face] = e2;
        precomp.vert0[face] = v0;
        precomp.pvec[face] = pv * i_det;
        precomp.is_parallel[face] = fabs(det) < epsilon;
    }
}

// Moller-Trumbore algorithm to find whether the ray and the triangle intersect or not
const bool Mesher::ray_intersects_triangle(const Vec3d &origin, const Vec3d &direction, const int face) {
    const double zero = -1.0 * epsilon; // Ray little bit outside the triangle is also OK
    const double one = 1.0 + epsilon;

    Vec3d s, q;
    double u, v;

    // Check if ray and triangle are effectively parallel
    // Must be in the beginning because the following calculations may become in-reliable (division by zero etc)
    if (precomp.is_parallel[face]) return false; // Comment out when using ray_intersects_surface(...)

    s = origin - precomp.vert0[face];
    u = s.dotProduct(precomp.pvec[face]);

    // Check if ray intersects triangle
    if (u < zero || u > one) return false;

    q = s.crossProduct(precomp.edge1[face]);
    v = direction.dotProduct(q); // if ray == [0,0,1] then v = q.z

    if (v < zero || u + v > one) return false;

    // Check whether there is a line intersection or ray intersection
    if (precomp.edge2[face].dotProduct(q) < 0) return false;

    return true;
}

// Determine how many times the ray intersects with triangles on the surface mesh
const int Mesher::ray_intersects_surface(Mesh* mesh, const Vec3d &origin, const Vec3d &direction) {
    int face, n0, n1, n2;
    int crossings = 0;

    vector<bool> node_passed(mesh->get_n_nodes(), false);

    // Loop through all the faces
    for (face = 0; face < mesh->get_n_faces(); ++face) {
        if (precomp.is_parallel[face]) continue;

        n0 = mesh->get_face(face, 0);
        n1 = mesh->get_face(face, 1);
        n2 = mesh->get_face(face, 2);

        // If the ray has already passed neighbouring triangle, skip current one
        if (node_passed[n0] || node_passed[n1] || node_passed[n2]) continue;

        if (ray_intersects_triangle(origin, direction, face)) {
            crossings++;
            if (crossings >= 2) return crossings;
            node_passed[n0] = true;
            node_passed[n1] = true;
            node_passed[n2] = true;
        }
    }

    return crossings;
}

// Determine whether the ray intersects with any of the triangles on the surface mesh
const bool Mesher::ray_intersects_surface_fast(Mesh* mesh, const Vec3d &origin, const Vec3d &direction) {
    // Loop through all the faces
    for (int face = 0; face < mesh->get_n_faces(); ++face) {
        if (ray_intersects_triangle(origin, direction, face))
            return true;
    }

    return false;
}

const void Mesher::mark_elems(Mesh* mesh, const AtomReader::Types* types) {
    int elem, node, location;
    int n_elems = mesh->get_n_elems();
    mesh->init_elemmarkers(n_elems);

    // Loop through all the elements
    for (elem = 0; elem < n_elems; ++elem) {
        location = 0;
        // Loop through all the nodes in element
        for (node = 0; node < n_nodes_per_elem; ++node) {
            int nodemarker = mesh->get_nodemarker(mesh->get_elem(elem, node));
            if (nodemarker == types->type_vacuum)
                location += 3; // If element has node both in vacuum and bulk, mark it as vacuum
            else if (nodemarker == types->type_bulk)
                location--;
        }
        // Element in vacuum is supposed not to have nodes in bulk
        if (location > 0)
            mesh->add_elemmarker(types->type_vacuum);
        // Element in bulk is supposed not to have nodes in vacuum
        else if (location < 0)
            mesh->add_elemmarker(types->type_bulk);
        // Element in surface is supposed consist only of nodes in surface
        else
            mesh->add_elemmarker(types->type_surf);
    }
}

/* Function to mark nodes with ray-triangle intersection technique
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle
 *
 * When the ray from the node in z-direction crosses the surface, the node is located in material,
 * otherwise it's in vacuum.
 * The technique works perfectly with surface without "mushroomish" tips. The mushrooms give
 * some false nodes but because of their low spatial density they can be ignored.
 */
const void Mesher::mark_nodes(Mesh* mesh, const AtomReader::Types* types, const bool postprocess) {
    int node, elem, coord;
    Vec3d ray_origin, ray_direction(0, 0, 1);

    mesh->init_nodemarkers(mesh->get_n_nodes());
    precompute_triangles(mesh, ray_direction);

    // Mark surface nodes by their known position in array
    for (node = mesh->indxs.surf_start; node <= mesh->indxs.surf_end; ++node)
        mesh->add_nodemarker(types->type_surf);

    // Mark bulk nodes by their known position in array
    for (node = mesh->indxs.bulk_start; node <= mesh->indxs.bulk_end; ++node)
        mesh->add_nodemarker(types->type_bulk);

    // Mark vacuum nodes by their known position in array
    for (node = mesh->indxs.vacuum_start; node <= mesh->indxs.vacuum_end; ++node)
        mesh->add_nodemarker(types->type_vacuum);

    // Loop through all the nodes made by tetgen
    // and mark them by vertical-ray-intersects-surface-faces-technique
    for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node) {
        // Compose the ray_origin vector
        for (coord = 0; coord < n_coordinates; ++coord)
            ray_origin[coord] = mesh->get_node(node, coord);

        // If ray intersects at least one surface triangle, the node
        // is considered to be located in bulk, otherwise it's located to vacuum
        if ( ray_intersects_surface_fast(mesh, ray_origin, ray_direction) )
            mesh->add_nodemarker(types->type_bulk);
        else
            mesh->add_nodemarker(types->type_vacuum);
    }

    // Mark the elements by the node markers
    mark_elems(mesh, types);

    // Post process the nodes in shadow areas
    if (postprocess)
        post_process_node_marking(mesh, types);
}

const void Mesher::post_process_node_marking(Mesh* mesh, const AtomReader::Types* types) {
    int elem, node, location;
    int n_elems = mesh->get_n_elems();

    bool node_changed = true;
    for (int safety_cntr = 0; safety_cntr < 10, node_changed; ++safety_cntr) {
        // Post-process the wrongly marked nodes in vacuum
        node_changed = false;
        for (elem = 0; elem < n_elems; ++elem) {
            if (mesh->get_elemmarker(elem) != types->type_vacuum) continue;

            // Force all the nodes in vacuum elements to be non-bulk ones
            for (node = 0; node < n_nodes_per_elem; ++node) {
                int n = mesh->get_elem(elem, node);
                if(mesh->get_nodemarker(n) == types->type_bulk) {
                    mesh->set_nodemarker(n, types->type_vacuum);
                    node_changed = true;
                }
            }
        }

        // If some of the nodes were changed,
        // mark all the surface and bulk elements again
        if(node_changed)
            for (elem = 0; elem < n_elems; ++elem) {
                if(mesh->get_elemmarker(elem) == types->type_vacuum) continue;

                location = 0;
                // Loop through all the nodes in element
                for (node = 0; node < n_nodes_per_elem; ++node) {
                    int nodemarker = mesh->get_nodemarker(mesh->get_elem(elem, node));
                    if (nodemarker == types->type_vacuum)
                        location += 3; // If element has node both in vacuum and bulk, mark it as vacuum
                    else if(nodemarker == types->type_bulk)
                        location--;
                }

                if (location > 0)
                    mesh->set_elemmarker(elem, types->type_vacuum);
                else if (location < 0)
                    mesh->set_elemmarker(elem, types->type_bulk);
                else
                    mesh->set_elemmarker(elem, types->type_surf);
            }
    }

    mesh->calc_statistics(types);
    require(mesh->stat.n_bulk > 0, "Nodemarker post processing deleted all the bulk atoms.\n"
            "Consider altering the surface refinment factor.");
}

/* Function to mark nodes with ray-triangle intersection technique
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle
 *
 * When the ray from the node in "random" direction crosses the surface 2*n times, the node
 * is located in vacuum, when it crosses the surface for 2*n+1, times, it's in material.
 * The technique works well with smooth surface. The non-smooth one gives for some ray direction
 * too many crossings. The workaround is to change the ray direction for the problematic nodes
 * to make sure the ray crosses only "quality" triangles.
 */
const void Mesher::mark_nodes_long(Mesh* mesh, const AtomReader::Types* types) {
    int node, coord;
    Vec3d ray_origin, ray_direction(0, 0, 1);
    vector<int> markers(mesh->get_n_nodes(), types->type_none);

    // Set of x & y directions for probing ray
    const double tilt_x[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!
    const double tilt_y[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!

    mesh->init_nodemarkers(mesh->get_n_nodes());

    // Mark surface nodes by their known position in array
    for (node = mesh->indxs.surf_start; node <= mesh->indxs.surf_end; ++node)
        mesh->add_nodemarker(types->type_surf);

    // Mark bulk nodes by their known position in array
    for (node = mesh->indxs.bulk_start; node <= mesh->indxs.bulk_end; ++node)
        mesh->add_nodemarker(types->type_bulk);

    // Mark vacuum nodes by their known position in array
    for (node = mesh->indxs.vacuum_start; node <= mesh->indxs.vacuum_end; ++node)
        mesh->add_nodemarker(types->type_vacuum);

    bool all_nodes_done = false;

    // Try different ray directions until all the nodes have <= 1 intersection with surface faces
    for (const int tx : tilt_x) {
        if (all_nodes_done) break;

        ray_direction[0] = tx;
        for (const int ty : tilt_x) {
            if (all_nodes_done) break;

            ray_direction[1] = ty;
            ray_direction.normalize();

            precompute_triangles(mesh, ray_direction);

            all_nodes_done = true;

            // Loop through all the nodes inserted by tetgen
            for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node) {
                // If node is already checked, skip it
                if (markers[node] != types->type_none) continue;

                // Compose the ray_origin vector
                for (coord = 0; coord < n_coordinates; ++coord)
                    ray_origin[coord] = mesh->get_node(node, coord);

                // Find how many times the ray from node crosses the surface
                int crossings = ray_intersects_surface_fast(mesh, ray_origin, ray_direction);

                if (crossings == 0)
                    markers[node] = types->type_vacuum;
                else if (crossings == 1)
                    markers[node] = types->type_bulk;
                else
                    all_nodes_done = false;
            }
        }
    }

    // Insert found markers to mesh
    for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node)
        mesh->add_nodemarker(markers[node]);
}

// Mark faces by the sequence of nodes
const void Mesher::mark_faces_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);
    bool difs[n_nodes_per_face];

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < n_nodes_per_face; ++j)
            difs[j] = mesh->get_face(i, j) < nmax;

        if (difs[0] + difs[1] + difs[2] == 3)
            mesh->add_facemarker(types->type_surf);
        else
            mesh->add_facemarker(types->type_none);
    }
}

// Function to determine whether the center of face is on the boundary of simulation cell or not
inline bool on_boundary(const double face, const double face_max) {
    const double eps = 0.1;
    return fabs(face - face_max) < eps;
}

// Mark the boundary faces of mesh
const void Mesher::mark_faces(Mesh* mesh, const AtomReader::Sizes* sizes,
        const AtomReader::Types* types) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);

    for (int i = 0; i < N; ++i) {
        if (on_boundary(mesh->get_face_centre(i, 0), sizes->xmin))
            mesh->add_facemarker(types->type_xmin);
        else if (on_boundary(mesh->get_face_centre(i, 0), sizes->xmax))
            mesh->add_facemarker(types->type_xmax);
        else if (on_boundary(mesh->get_face_centre(i, 1), sizes->ymin))
            mesh->add_facemarker(types->type_ymin);
        else if (on_boundary(mesh->get_face_centre(i, 1), sizes->ymax))
            mesh->add_facemarker(types->type_ymax);
        else if (on_boundary(mesh->get_face_centre(i, 2), sizes->zminbox))
            mesh->add_facemarker(types->type_zmin);
        else if (on_boundary(mesh->get_face_centre(i, 2), sizes->zmaxbox))
            mesh->add_facemarker(types->type_zmax);
        else
            mesh->add_facemarker(types->type_surf);
    }
}

// Function to remove too big tetrahedra from the mesh
const void Mesher::clean_elems(Mesh* mesh, const double rmax, const string cmd) {
    const REAL* node = mesh->get_nodes();      // pointer to nodes
    const int* elem = mesh->get_elems(); // pointer to tetrahedron corners

    double dx, r12, r13, r14, r23, r24, r34;
    int i, j, k, l1, l2, l3, l4;

    int N = mesh->get_n_elems();
    vector<bool> isQuality(N);

    // Loop through the tetrahedra
    for (i = 0; i < N; ++i) {
        j = 4 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r14 = r23 = r24 = r34 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * elem[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * elem[j + 1] + k;   // ..2nd node
            l3 = 3 * elem[j + 2] + k; // ..3rd node
            l4 = 3 * elem[j + 3] + k; // ..4th node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of tetrahedron 1st edge
            dx = node[l1] - node[l3];
            r13 += dx * dx; // ..2nd edge
            dx = node[l1] - node[l4];
            r14 += dx * dx; // ..3rd edge
            dx = node[l2] - node[l3];
            r23 += dx * dx; // ..4th edge
            dx = node[l2] - node[l4];
            r24 += dx * dx; // ..5th edge
            dx = node[l3] - node[l4];
            r34 += dx * dx; // ..6th edge
        }
        // Keep only the elements with appropriate quality
        isQuality[i] = (r12 < rmax) & (r13 < rmax) & (r14 < rmax) & (r23 < rmax) & (r24 < rmax)
                & (r34 < rmax);
    }

    // Make temporary copy from the input mesh
    shared_ptr<Mesh> temp_mesh(new Mesh(this->mesher));
    temp_mesh->init_elems(mesh->get_n_elems());
    temp_mesh->copy_elems(mesh, 0);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_elems(vector_sum(isQuality));

    // Insert tetrahedra with suitable quality
    update_list((int*)mesh->get_elems(), temp_mesh->get_elems(), isQuality, 4);

    mesh->recalc(cmd);
}

// Function to remove too big faces from the mesh
const void Mesher::clean_faces(Mesh* mesh, const double rmax, const string cmd) {
    const REAL* node = mesh->get_nodes();		// pointer to nodes
    const int* face = mesh->get_faces();		// pointer to triangular faces

    double dx, r12, r13, r23;
    int i, j, k, l1, l2, l3;
    int N = mesh->get_n_faces();
    vector<bool> isQuality(N);

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r23 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k;	// ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of face's 1st edge
            dx = node[l1] - node[l3];
            r13 += dx * dx; // ..2nd edge
            dx = node[l2] - node[l3];
            r23 += dx * dx; // ..3rd edge
        }
        // Keep only the faces with appropriate quality
        isQuality[i] = (r12 < rmax) & (r13 < rmax) & (r23 < rmax);
    }

    // Make temporary copy from the input mesh
    shared_ptr<Mesh> temp_mesh(new Mesh(this->mesher));
    temp_mesh->init_faces(mesh->get_n_faces());
    temp_mesh->copy_faces(mesh, 0);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_faces(vector_sum(isQuality));

    // Insert faces with suitable quality
    update_list((int*)mesh->get_faces(), temp_mesh->get_faces(), isQuality, 3);

    mesh->recalc(cmd);
}

// Function to remove the objects from element, face or edge list
const void Mesher::update_list(int* new_list, const int* old_list, const vector<bool> is_quality, int M) {
    int i, j, k;
    int N = is_quality.size();
    j = 0;
    // loop through the old array and move quality objects to the left
    for (i = 0; i < N; ++i)
        if (is_quality[i]) {
            for (k = 0; k < M; ++k)
                new_list[j + k] = old_list[M * i + k];
            j += M;
        }
}

} /* namespace femocs */
