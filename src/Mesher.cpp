/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

#include <memory>
#include <algorithm>
#include <time.h>
#include <stdlib.h>     /* srand, rand */

using namespace std;
namespace femocs {

Mesher::Mesher(const string mesher, const double latconst) {
    this->latconst = latconst;
    this->mesher = mesher;
}

// Function to generate simple mesh that consists of one tetrahedron
const void Mesher::get_test_mesh(Mesh* new_mesh) {
    const int n_nodes = 4;
    const int n_faces = 4;
    const int n_elems = 1;

    new_mesh->init_nodes(n_nodes);
    new_mesh->add_node(1.0, 0.0, 0.7);
    new_mesh->add_node(-1.0, 0.0, 0.7);
    new_mesh->add_node(0.0, 1.0, -0.7);
    new_mesh->add_node(0.0, -1.0, -0.7);

    new_mesh->init_faces(n_faces);
    new_mesh->add_face(0, 1, 3);
    new_mesh->add_face(1, 2, 3);
    new_mesh->add_face(2, 0, 3);
    new_mesh->add_face(0, 1, 2);

    new_mesh->init_elems(n_elems);
    new_mesh->add_elem(0, 1, 2, 3);
}

// Function to generate mesh from bulk and vacuum atoms
const void Mesher::get_volume_mesh(Mesh* new_mesh, Bulk* bulk, Surface* surf, Vacuum* vacuum,
        const string cmd) {
    int i;
    int n_bulk = bulk->get_n_atoms();
    int n_surf = surf->get_n_atoms();
    int n_vacuum = vacuum->get_n_atoms();

    new_mesh->init_nodes(n_bulk + n_surf + n_vacuum);

    for (i = 0; i < n_surf; ++i)
        new_mesh->add_node(surf->get_x(i), surf->get_y(i), surf->get_z(i));

    for (i = 0; i < n_bulk; ++i)
        new_mesh->add_node(bulk->get_x(i), bulk->get_y(i), bulk->get_z(i));

    for (i = 0; i < n_vacuum; ++i)
        new_mesh->add_node(vacuum->get_x(i), vacuum->get_y(i), vacuum->get_z(i));

    new_mesh->recalc(cmd);
}

const void Mesher::extract_mesh(vector<bool>* is_vacuum, Mesh* big_mesh, const int n_bulk,
        const int n_surf, const double zmin) {
    int i, j;
    int n_elems = big_mesh->get_n_elems();

    is_vacuum->resize(n_elems);

    int n_vacuum, n_not_bottom, n_not_bulk;

    for (i = 0; i < n_elems; ++i) {
        n_vacuum = n_not_bottom = n_not_bulk = 0;

        // Find nodes that are not inside the bulk
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (big_mesh->get_elem(i, j) >= (n_bulk + n_surf)) n_vacuum++;

        // Find nodes that are not in floating element
        for (j = 0; j < n_nodes_per_elem; ++j)
            if ((big_mesh->get_elem(i, j) >= n_bulk) || (big_mesh->get_elem(i, j) < n_surf))
                n_not_bulk++;

        // Find elements not belonging to the bottom of simulation cell
        for (j = 0; j < n_nodes_per_elem; ++j) {
            int node = big_mesh->get_elem(i, j);
            if (fabs(big_mesh->get_z(node) - zmin) > 1.0 * latconst) n_not_bottom++;
        }

        //(*is_vacuum)[i] = (n_vacuum > 0) && (n_not_bulk > 1) && (n_not_bottom == n_nodes_per_elem);
        (*is_vacuum)[i] = (n_vacuum > 0) && (n_not_bulk > 1);
    }
}

const void Mesher::extract_mesh_bymarker(vector<bool>* is_vacuum, Mesh* big_mesh, const AtomReader::Types* types) {
    int elem, node;
    int n_elems = big_mesh->get_n_elems();

    is_vacuum->resize(n_elems);

    int n_vacuum, n_not_bulk;

    for (elem = 0; elem < n_elems; ++elem) {
        n_vacuum = n_not_bulk = 0;

        // Find nodes that are not inside the bulk
        for (node = 0; node < n_nodes_per_elem; ++node)
            if ( big_mesh->get_nodemarker(big_mesh->get_elem(elem,node)) == types->type_vacuum )
                n_vacuum++;

        // Find nodes that are in the surface
        for (node = 0; node < n_nodes_per_elem; ++node)
            if ( big_mesh->get_nodemarker(big_mesh->get_elem(elem,node)) != types->type_bulk )
                n_not_bulk++;

        (*is_vacuum)[elem] = (n_vacuum > 0) && (n_not_bulk > 1);
    }
}

const void Mesher::separate_vacuum_mesh(Mesh* vacuum_mesh, Mesh* big_mesh, const int n_bulk,
        const int n_surf, const double zmin, const string cmd) {
    int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, n_bulk, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    vacuum_mesh->init_nodes(big_mesh->get_n_nodes());
    vacuum_mesh->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    vacuum_mesh->init_elems(n_vacuum_elems);

    for (int i = 0; i < big_mesh->get_n_elems(); ++i)
        if (elem_in_vacuum[i]) {
            int j = n_nodes_per_elem * i;
            vacuum_mesh->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    vacuum_mesh->recalc(cmd);
}

const void Mesher::separate_bulk_mesh(Mesh* bulk, Mesh* big_mesh, const int n_bulk,
        const int n_surf, const double zmin, const string cmd) {
    int* elems = big_mesh->get_elems();
    int n_elems = big_mesh->get_n_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, n_bulk, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    bulk->init_nodes(big_mesh->get_n_nodes());
    bulk->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    bulk->init_elems(n_elems - n_vacuum_elems);

    for (int i = 0; i < n_elems; ++i)
        if (!elem_in_vacuum[i]) {
            int j = n_nodes_per_elem * i;
            bulk->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    bulk->recalc(cmd);
}

const void Mesher::separate_meshes(Mesh* bulk, Mesh* vacuum, Mesh* big_mesh, const int n_bulk,
        const int n_surf, const double zmin, const string cmd) {

    int i, j;
    int n_nodes = big_mesh->get_n_nodes();
    int n_elems = big_mesh->get_n_elems();
    int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, n_bulk, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

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
    int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh_bymarker(&elem_in_vacuum, big_mesh, types);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

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
const void Mesher::generate_monolayer_surf_faces(Mesh* mesh, const int n_bulk, const int n_surf) {
    int i, j, n1, n2, n3;
    int locs[n_nodes_per_elem];
    bool b_locs[n_nodes_per_elem];

    int n_elems = mesh->get_n_elems();
    int n_nodes = mesh->get_n_nodes();
    int n_faces = mesh->get_n_faces();

    int* elems = mesh->get_elems();
    vector<bool> is_surface(n_elems);

    // Mark the elements that have exactly one face on the surface and one node in the vacuum
    for (i = 0; i < n_elems; ++i) {
        for (j = 0; j < n_nodes_per_elem; ++j) {
            int elem = mesh->get_elem(i, j);
            // Node in bulk?
            if ((elem >= n_surf) && (elem < (n_surf + n_bulk)))
                locs[j] = -1;

            // Node in vacuum?
            else if (elem > (n_bulk + n_surf))
                locs[j] = 0;

            // Node in surface!
            else
                locs[j] = 1;
        }

        // Bulk, surf and vacuum are encoded in such a way that
        // sum(locs) == 3 only if 3 nodes are on surface and 1 in vacuum
        is_surface[i] = locs[0] + locs[1] + locs[2] + locs[3] == 3;
    }

    int n_surf_faces = accumulate(is_surface.begin(), is_surface.end(), 0);

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
                b_locs[j] = mesh->get_elem(i, j) <= (n_bulk + n_surf);

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
const void Mesher::generate_surf_faces(Mesh* mesh, const int n_surf) {
    int elem, node, n0, n1, n2;

    int n_elems = mesh->get_n_elems();
    int* elems = mesh->get_elems();

    vector<bool> elem_in_surface(n_elems);
    bool surf_locs[n_nodes_per_elem];

    // Mark the elements that have exactly one face on the surface and one node in the vacuum
    for (elem = 0; elem < n_elems; ++elem) {
        for (node = 0; node < n_nodes_per_elem; ++node)
            surf_locs[node] = mesh->get_elem(elem, node) < n_surf;

        elem_in_surface[elem] = surf_locs[0] + surf_locs[1] + surf_locs[2] + surf_locs[3] == 3;
    }

    int n_surf_faces = accumulate(elem_in_surface.begin(), elem_in_surface.end(), 0);
    mesh->init_faces(n_surf_faces);

    // Generate the faces that separate material and vacuum
    // It is assumed that each element has only one face on the surface
    for (elem = 0; elem < n_elems; ++elem)
        if (elem_in_surface[elem]) {
            for (node = 0; node < n_nodes_per_elem; ++node)
                surf_locs[node] = mesh->get_elem(elem, node) < n_surf;

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

void Mesher::precompute_triangles(Mesh* mesh, const Vec3d &direction) {
    Vec3d v0, v1, v2, e1, e2, pv;
    int face, coord;
    double det, i_det;
    int n_faces = mesh->get_n_faces();

    edge1.reserve(n_faces);
    edge2.reserve(n_faces);
    vert0.reserve(n_faces);
    pvec.reserve(n_faces);
    is_parallel.reserve(n_faces);

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

        edge1[face] = e1 * i_det;
        edge2[face] = e2;
        vert0[face] = v0;
        pvec[face] = pv * i_det;
        is_parallel[face] = fabs(det) < epsilon;
    }
}

// Moller-Trumbore algorithm to find whether the ray and the triangle intersect or not
bool Mesher::ray_intersects_triangle(const Vec3d &origin, const Vec3d &direction, const int face) {
    const double zero = -1.0 * epsilon; // Ray little bit outside the triangle is also OK
    const double one = 1.0 + epsilon;

    Vec3d s, q;
    double u, v;

    if (is_parallel[face]) return false; // Comment in when using ray_surface_intersects_fast

    s = origin - vert0[face];
    u = s.dotProduct(pvec[face]);

    if (u < zero || u > one) return false;

    q = s.crossProduct(edge1[face]);
    v = direction.dotProduct(q); // if ray == [0,0,1] then v = q.z

    if (v < zero || u + v > one) return false;

    // Check whether there is a line intersection or ray intersection
    if (edge2[face].dotProduct(q) < 0) return false;

    return true;
}

// Determine how many times the ray intersects with triangles on the surface mesh
int Mesher::ray_surface_intersects(Mesh* mesh, const Vec3d &origin, const Vec3d &direction) {
    int face, n0, n1, n2;
    int crossings = 0;

    vector<bool> node_passed(mesh->get_n_nodes(), false);

    // Loop through all the faces
    for (face = 0; face < mesh->get_n_faces(); ++face) {
        if (is_parallel[face]) continue;

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
bool Mesher::ray_surface_intersects_fast(Mesh* mesh, const Vec3d &origin, const Vec3d &direction) {
    // Loop through all the faces
    for (int face = 0; face < mesh->get_n_faces(); ++face) {
        if (ray_intersects_triangle(origin, direction, face))
            return true;
    }

    return false;
}

/* Function to mark nodes with ray-triangle intersection technique
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle
 *
 * When the ray from the node in z-direction crosses the surface, the node is located in material,
 * otherwise it's in vacuum.
 * The technique works perfectly with surface without "mushroomish" tips. The mushrooms give
 * some false nodes but because of their low spatial density they can be ignored.
 */
void Mesher::mark_nodes(Mesh* mesh, const AtomReader::Types* types, const int n_surf) {
    int node, coord;
    Vec3d ray_origin, ray_direction(0, 0, 1);

    mesh->init_nodemarkers(mesh->get_n_nodes());

    // Mark surface nodes by their position in array
    for (node = 0; node < n_surf; ++node)
        mesh->add_nodemarker(types->type_surf);

    precompute_triangles(mesh, ray_direction);

    // Loop through all the non-surface nodes
    for (node = n_surf; node < mesh->get_n_nodes(); ++node) {
        // Compose the ray_origin vector
        for (coord = 0; coord < n_coordinates; ++coord)
            ray_origin[coord] = mesh->get_node(node, coord);

        if ( ray_surface_intersects_fast(mesh, ray_origin, ray_direction) )
            mesh->add_nodemarker(types->type_bulk);
        else
            mesh->add_nodemarker(types->type_vacuum);
    }
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
void Mesher::mark_nodes_long(Mesh* mesh, const AtomReader::Types* types, const int n_surf) {
    int node, coord;
    Vec3d ray_origin, ray_direction(0, 0, 1);
    vector<int> markers(mesh->get_n_nodes(), types->type_none);

    // Set of x & y directions for probing ray
    const double tilt_x[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!
    const double tilt_y[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!

    mesh->init_nodemarkers(mesh->get_n_nodes());

    // Mark surface nodes by their position in array
    for (node = 0; node < n_surf; ++node)
        mesh->add_nodemarker(types->type_surf);

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

            // Loop through all the non-surface nodes
            for (node = n_surf; node < mesh->get_n_nodes(); ++node) {
                // If node is already checked, skip it
                if (markers[node] != types->type_none) continue;

                // Compose the ray_origin vector
                for (coord = 0; coord < n_coordinates; ++coord)
                    ray_origin[coord] = mesh->get_node(node, coord);

                // Find how many times the ray from node crosses the surface
                int crossings = ray_surface_intersects_fast(mesh, ray_origin, ray_direction);

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
    for (node = n_surf; node < mesh->get_n_nodes(); ++node)
        mesh->add_nodemarker(markers[node]);
}

/* Mark elements by the sequence of nodes
 * The idea is that when Tetgen adds nodes to the mesh, it adds them to the end of node list.
 * In that ways the first nodes will always be the surface and bulk nodes and there is no need
 * to search them later on.
 */
void Mesher::mark_elems_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types) {
    int N = mesh->get_n_elems();
    mesh->init_elemmarkers(N);
    bool difs[4];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 4; ++j)
            difs[j] = mesh->get_elem(i, j) < nmax;

        if (difs[0] + difs[1] + difs[2] + difs[3] == 3) {
            if (difs[0] & difs[3])
                mesh->add_elemmarker(1);
            else if (difs[0] & (!difs[3]))
                mesh->add_elemmarker(2);
            else if ((!difs[0]) & difs[3])
                mesh->add_elemmarker(2);
            else if ((!difs[0]) & (!difs[3])) mesh->add_elemmarker(4);
        } else
            mesh->add_elemmarker(types->type_none);
    }
}

// Mark faces by the sequence of nodes
void Mesher::mark_faces_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types) {
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
void Mesher::mark_faces(Mesh* mesh, const AtomReader::Sizes* sizes,
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
void Mesher::clean_elems(Mesh* mesh, const double rmax, const string cmd) {
    REAL* node = mesh->get_nodes();      // pointer to nodes
    int* elem = mesh->get_elems(); // pointer to tetrahedron corners

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
    mesh->init_elems(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert tetrahedra with suitable quality
    update_list(mesh->get_elems(), temp_mesh->get_elems(), isQuality, 4);

    mesh->recalc(cmd);
}

// Function to remove too big faces from the mesh
void Mesher::clean_faces(Mesh* mesh, const double rmax, const string cmd) {
    REAL* node = mesh->get_nodes();		// pointer to nodes
    int* face = mesh->get_faces();		// pointer to triangular faces

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
    mesh->init_faces(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert faces with suitable quality
    update_list(mesh->get_faces(), temp_mesh->get_faces(), isQuality, 3);

    mesh->recalc(cmd);
}

// Function to remove the objects from element, face or edge list
void Mesher::update_list(int* new_list, int* old_list, vector<bool> is_quality, int M) {
    int i, j, k;
    int N = is_quality.size();
    j = 0;
    // loop through the old array and move quality objects to the left
    for (i = 0; i < N; ++i) {
        if (is_quality[i]) {
            for (k = 0; k < M; ++k)
                new_list[j + k] = old_list[M * i + k];
            j += M;
        }
    }
}

} /* namespace femocs */
