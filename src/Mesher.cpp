/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"
#include <fstream>

using namespace std;
namespace femocs {

/* =======================================
 *  Implementation of RaySurfaceIntersect
 * ======================================= */

RaySurfaceIntersect::RaySurfaceIntersect(Mesh* mesh) {
    reserve(0);
    faces.reserve(0);
    this->mesh = mesh;
}

// Reserve memory for precompute data
const void RaySurfaceIntersect::reserve(const int n) {
    edge1.reserve(n);
    edge2.reserve(n);
    vert0.reserve(n);
    pvec.reserve(n);
    is_parallel.reserve(n);
}

// Precompute the data needed to execute the Moller-Trumbore algorithm
const void RaySurfaceIntersect::precompute_triangles(const Vec3 &direction) {
    const int n_faces = faces.size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        SimpleFace sface = faces[face];
        Vec3 v0 = mesh->get_vec(sface[0]);
        Vec3 v1 = mesh->get_vec(sface[1]);
        Vec3 v2 = mesh->get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = direction.crossProduct(e2);
        double det = e1.dotProduct(pv);
        double i_det = 1.0 / det;

        edge1[face] = e1 * i_det;
        edge2[face] = e2;
        vert0[face] = v0;
        pvec[face] = pv * i_det;
        is_parallel[face] = fabs(det) < epsilon;
    }

    // Calculate the min & max of x, y & z coordinates
    mesh->calc_statistics();
    // Store the coordinates of surface corner node
    corner_node = mesh->get_node(0);
}

// Moller-Trumbore algorithm to find whether the ray and the triangle intersect or not
const bool RaySurfaceIntersect::ray_intersects_triangle(const Vec3 &origin, const Vec3 &direction,
        const int face) {
    Vec3 s, q;
    double u, v;

    // Check if ray and triangle are effectively parallel
    // Must be in the beginning because the following calculations may become in-reliable (division by zero etc)
    if (is_parallel[face]) return false; // Comment out when using ray_intersects_surface(...)

    s = origin - vert0[face];
    u = s.dotProduct(pvec[face]);

    // Check if ray intersects triangle
    // Ray little bit outside the triangle is also OK
    if (u < zero || u > one) return false;

    q = s.crossProduct(edge1[face]);
    v = direction.dotProduct(q); // if ray == [0,0,1] then v = q.z

    if (v < zero || u + v > one) return false;

    // Check whether there is a line intersection or ray intersection
    if (edge2[face].dotProduct(q) < 0) return false;

    return true;
}

// Determine how many times the ray intersects with triangles on the surface mesh
const int RaySurfaceIntersect::ray_intersects_surface(const Vec3 &origin, const Vec3 &direction) {
    const int n_faces = faces.size();
    int crossings = 0;

    vector<bool> node_passed(mesh->get_n_nodes(), false);

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        if (is_parallel[face]) continue;

        SimpleFace sface = faces[face];

        // If the ray has already passed neighbouring triangle, skip current one
        if (node_passed[sface[0]] || node_passed[sface[1]] || node_passed[sface[2]]) continue;

        if (ray_intersects_triangle(origin, direction, face)) {
            crossings++;
            if (crossings >= 2) return crossings;
            node_passed[sface[0]] = true;
            node_passed[sface[1]] = true;
            node_passed[sface[2]] = true;
        }
    }

    return crossings;
}

// Determine whether the ray intersects with any of the triangles on the surface mesh
const bool RaySurfaceIntersect::ray_intersects_surface_fast(const Vec3 &origin,
        const Vec3 &direction) {
    const int n_faces = faces.size();

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        if (ray_intersects_triangle(origin, direction, face)) return true;
    }

    return false;
}

// Function to manually generate surface faces from elements and surface nodes
const void RaySurfaceIntersect::generate_surf_faces() {
    int elem;

    const int n_elems = mesh->get_n_elems();
    const int max_surf_indx = mesh->indxs.surf_end;
    const int* elems = mesh->get_elems();

    // booleans showing whether element i has exactly one face on the surface or not
    vector<bool> elem_in_surface(n_elems);
    // booleans showing whether node i is on the surface or not
    vector<bool> surf_locs;

    // Mark the elements that have exactly one face on the surface
    for (elem = 0; elem < n_elems; ++elem) {
        surf_locs = mesh->get_simpleelem(elem) <= max_surf_indx;
        elem_in_surface[elem] = (vector_sum(surf_locs) == 3);
    }

    // Reserve memory for surface faces
    faces.reserve(2 + vector_sum(elem_in_surface));

    // Make two big faces that pass the xy-min-max corners of surface
    faces.push_back(SimpleFace(0, 1, 2));
    faces.push_back(SimpleFace(0, 2, 3));

    // Generate the faces that separate material and vacuum
    // The faces are taken from the elements that have exactly one face on the surface
    for (elem = 0; elem < n_elems; ++elem)
        if (elem_in_surface[elem]) {
            // Find the indices of nodes that are on the surface
            surf_locs = mesh->get_simpleelem(elem) <= max_surf_indx;

            /* The possible combinations of surf_locs and n0,n2,n3:
             * surf_locs: 1110   1101   1011   0111
             *        n0: elems0 elems0 elems0 elems1
             *        n1: elems1 elems1 elems2 elems2
             *        n2: elems2 elems3 elems3 elems3
             */
            int node = mesh->n_nodes_per_elem * elem;
            int n0 = surf_locs[0] * elems[node + 0] + (!surf_locs[0]) * elems[node + 1];
            int n1 = (surf_locs[0] & surf_locs[1]) * elems[node + 1]
                    + (surf_locs[2] & surf_locs[3]) * elems[node + 2];
            int n2 = (!surf_locs[3]) * elems[node + 2] + surf_locs[3] * elems[node + 3];
            faces.push_back(SimpleFace(n0, n1, n2));
        }
}

// Function to manually generate surface faces from elements, surface nodes and vacuum nodes
const void RaySurfaceIntersect::generate_monolayer_surf_faces() {
    int i, j;
    const int n_elems = mesh->get_n_elems();
    const int* elems = mesh->get_elems();
    vector<bool> is_surface(n_elems);

    // Mark the elements that have exactly one face on the surface and one node in the vacuum
    for (i = 0; i < n_elems; ++i) {

        int location = 0;
        for (int node : mesh->get_simpleelem(i)) {
            // Node in bulk?
            if ((node >= mesh->indxs.bulk_start) && (node <= mesh->indxs.bulk_end))
                location--;

            // Node in surface?
            else if (node <= mesh->indxs.surf_end)
                location++;

            // Node in vacuum!
//            else
//                location += 0;
        }

        // Bulk, surface and vacuum are encoded in such a way that
        // location == 3 only if 3 nodes are on surface and 1 in vacuum
        is_surface[i] = location == 3;
    }

    faces.reserve(vector_sum(is_surface));

    // Generate the faces that separate material and vacuum
    // It is assumed that each element has only one face on the surface
    for (i = 0; i < n_elems; ++i)
        if (is_surface[i]) {
            vector<bool> surf_locs = mesh->get_simpleelem(i) <= mesh->indxs.bulk_end;

            /* The possible combinations of b_locs and n1,n2,n3:
             * b_locs: 1110  1101  1011  0111
             *     n1: elem0 elem0 elem0 elem1
             *     n2: elem1 elem1 elem2 elem2
             *     n3: elem2 elem3 elem3 elem3
             */
            j = mesh->n_nodes_per_elem * i;
            int n1 = surf_locs[0] * elems[j + 0] + (!surf_locs[0]) * elems[j + 1];
            int n2 = (surf_locs[0] & surf_locs[1]) * elems[j + 1]
                    + (surf_locs[2] & surf_locs[3]) * elems[j + 2];
            int n3 = (!surf_locs[3]) * elems[j + 2] + surf_locs[3] * elems[j + 3];
            faces.push_back(SimpleFace(n1, n2, n3));
        }
}

// Function to output surface faces in .vtk format
const void RaySurfaceIntersect::write_faces(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    const string ftype = get_file_type(file_name);
    require(ftype == "vtk", "Unimplemented file type!");

    ofstream out(file_name);
    require(out.is_open(), "Can't open a file " + file_name);

    out.precision(8);

    out << "# vtk DataFile Version 3.0\n";
    out << "# Surface faces generated by RaySurfaceIntersect\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    const size_t n_nodes = mesh->get_n_nodes();
    const size_t n_faces = faces.size();
    const size_t celltype = 5; // 1-vertex, 5-triangle, 10-tetrahedron

    out << "POINTS " << n_nodes << " double\n";
    for (size_t node = 0; node < n_nodes; ++node)
        out << mesh->get_node(node) << "\n";

    out << "CELLS " << n_faces << " " << n_faces * (mesh->n_nodes_per_face + 1) << "\n";
    for (size_t face = 0; face < n_faces; ++face)
        out << mesh->n_nodes_per_face << " " << faces[face] << "\n";

    out << "CELL_TYPES " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; ++i)
        out << celltype << "\n";
}

/* ==========================
 *  Implementation of Mesher
 * ========================== */

Mesher::Mesher(Mesh* mesh) {
    this->mesh = mesh;
}

// Function to generate mesh from surface, bulk and vacuum atoms
const void Mesher::generate_mesh(Bulk &bulk, Surface &surf, Vacuum &vacuum, const string& cmd) {
    int i;
    int n_bulk = bulk.get_n_atoms();
    int n_surf = surf.get_n_atoms();
    int n_vacuum = vacuum.get_n_atoms();

    mesh->init_nodes(n_bulk + n_surf + n_vacuum);

    // Add surface atoms first,...
    for (i = 0; i < n_surf; ++i)
        mesh->add_node(surf.get_point(i));

    // ... bulk atoms second,...
    for (i = 0; i < n_bulk; ++i)
        mesh->add_node(bulk.get_point(i));

    // ... and vacuum atoms last
    for (i = 0; i < n_vacuum; ++i)
        mesh->add_node(vacuum.get_point(i));

    mesh->indxs.surf_start = 0;
    mesh->indxs.surf_end = mesh->indxs.surf_start + n_surf - 1;
    mesh->indxs.bulk_start = mesh->indxs.surf_end + 1;
    mesh->indxs.bulk_end = mesh->indxs.bulk_start + n_bulk - 1;
    mesh->indxs.vacuum_start = mesh->indxs.bulk_end + 1;
    mesh->indxs.vacuum_end = mesh->indxs.vacuum_start + n_vacuum - 1;
    mesh->indxs.tetgen_start = mesh->indxs.vacuum_end + 1;

    mesh->recalc("Q", cmd);
}

// Separate vacuum and bulk mesh from the union mesh by the node and element markers
const void Mesher::separate_meshes(Mesh* vacuum, Mesh* bulk, const string& cmd) {
    int i, j;
    const int n_nodes = mesh->get_n_nodes();
    const int n_elems = mesh->get_n_elems();

    const vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);
    const vector<bool> node_in_vacuum = vector_not(mesh->get_nodemarkers(), TYPES.BULK);
    const vector<bool> node_in_bulk = vector_not(mesh->get_nodemarkers(), TYPES.VACUUM);

    // Copy the non-bulk nodes from input mesh
    vacuum->copy_nodes(this->mesh, node_in_vacuum);
    // Generate mapping between old and new node indices
    vector<int> vacm_map(n_nodes, -1);
    for (i = 0, j = 0; i < n_nodes; ++i)
        if (node_in_vacuum[i]) vacm_map[i] = j++;

    // Copy the non-vacuum nodes from input mesh
    bulk->copy_nodes(this->mesh, node_in_bulk);
    // Generate mapping between old and new node indices
    vector<int> bulk_map(n_nodes, -1);
    for (i = 0, j = 0; i < n_nodes; ++i)
        if (node_in_bulk[i]) bulk_map[i] = j++;

    // Reserve memory for elements
    const int n_vacuum_elems = vector_sum(elem_in_vacuum);
    const int n_bulk_elems = n_elems - n_vacuum_elems;
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_bulk_elems);

    // Separate vacuum and bulk elements and arrange the node indices
    for (i = 0; i < n_elems; ++i) {
        SimpleElement se = mesh->get_simpleelem(i);
        if (elem_in_vacuum[i])
            vacuum->add_elem(vacm_map[se.n1], vacm_map[se.n2], vacm_map[se.n3], vacm_map[se.n4]);
        else
            bulk->add_elem(bulk_map[se.n1], bulk_map[se.n2], bulk_map[se.n3], bulk_map[se.n4]);
    }

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

// Separate vacuum mesh from the union mesh by the node and element markers
const void Mesher::separate_meshes(Mesh* vacuum, const string& cmd) {
    int i, j;
    const int n_nodes = mesh->get_n_nodes();
    const int n_elems = mesh->get_n_elems();

    const vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);
    const vector<bool> node_in_vacuum = vector_not(mesh->get_nodemarkers(), TYPES.BULK);

    // Copy the non-bulk nodes from input mesh
    vacuum->copy_nodes(this->mesh, node_in_vacuum);
    // Generate mapping between old and new node indices
    vector<int> vacm_map(n_nodes, -1);
    for (i = 0, j = 0; i < n_nodes; ++i)
        if (node_in_vacuum[i]) vacm_map[i] = j++;

    // Reserve memory for elements
    vacuum->init_elems(vector_sum(elem_in_vacuum));

    // Separate vacuum elements and arrange the node indices
    for (i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i]) {
            SimpleElement se = mesh->get_simpleelem(i);
            vacuum->add_elem(vacm_map[se.n1], vacm_map[se.n2], vacm_map[se.n3], vacm_map[se.n4]);
        }

    vacuum->recalc(cmd);
}

// Separate bulk and vacuum mesh from the union mesh by the element markers
const void Mesher::separate_meshes_noclean(Mesh* vacuum, Mesh* bulk, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    // Copy the non-vacuum nodes from input mesh without modification
    bulk->copy_nodes(this->mesh);

    // Reserve memory for elements
    const int n_vacuum_elems = vector_sum(elem_in_vacuum);
    const int n_bulk_elems = n_elems - n_vacuum_elems;
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_bulk_elems);

    // Separate vacuum and bulk elements
    for (int i = 0; i < n_elems; ++i) {
        if (elem_in_vacuum[i])
            vacuum->add_elem(mesh->get_simpleelem(i));
        else
            bulk->add_elem(mesh->get_simpleelem(i));
    }

    vacuum->recalc();
    bulk->recalc();
}

// Separate vacuum mesh from the union mesh by the element markers
const void Mesher::separate_meshes_noclean(Mesh* vacuum, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    // Reserve memory for elements
    vacuum->init_elems(vector_sum(elem_in_vacuum));

    // Separate vacuum elements
    for (int i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i])
            vacuum->add_elem(mesh->get_simpleelem(i));

    vacuum->recalc(cmd);
}

// Separate vacuum mesh from the union mesh by the element markers
// Handle nodes that belong only to one element
const void Mesher::separate_meshes_vol2(Mesh* vacuum, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const int n_nodes = mesh->get_n_nodes();
    vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);

    // Get the number of elements the node is connected with
    vector<int> n_elems_for_node(n_nodes, 0);

    for (int elem = 0; elem < n_elems; ++elem)
        if (elem_in_vacuum[elem]) {
            for (int node : mesh->get_simpleelem(elem))
                n_elems_for_node[node]++;
        }

    // Assign elements that are sharp in vacuum to bulk where they are blunt
   for (int elem = 0; elem < n_elems; ++elem)
       if (elem_in_vacuum[elem]) {
           bool elem_is_blunt = true;

           // element is sharp if at least one of its nodes belongs to only one element
           // nodes 0-3 are the corners of surface, so skip them because their sharpness is OK
           for (int node : mesh->get_simpleelem(elem))
               elem_is_blunt &= (n_elems_for_node[node] > 1) || (node < 4);

           elem_in_vacuum[elem] = elem_is_blunt ;
       }

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    // Reserve memory for elements
    vacuum->init_elems(vector_sum(elem_in_vacuum));

    // Separate vacuum elements
    for (int i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i])
            vacuum->add_elem(mesh->get_simpleelem(i));

    vacuum->recalc(cmd);
}

// Separate vacuum and bulk mesh from the union mesh by the element markers
// Handle nodes that belong only to one element
const void Mesher::separate_meshes_vol2(Mesh* vacuum, Mesh* bulk, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const int n_nodes = mesh->get_n_nodes();
    vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);

    // Get the number of elements the node is connected with
    vector<int> n_elems_for_node(n_nodes, 0);

    for (int elem = 0; elem < n_elems; ++elem)
        if (elem_in_vacuum[elem])
            for (int node : mesh->get_simpleelem(elem))
                n_elems_for_node[node]++;

    // Assign elements that are sharp in vacuum to bulk where they are blunt
   for (int elem = 0; elem < n_elems; ++elem)
       if (elem_in_vacuum[elem]) {
           bool elem_is_blunt = true;

           // element is sharp if at least one of its nodes belongs to only one element
           // nodes 0-3 are the corners of surface, so skip them because their sharpness is OK
           for (int node : mesh->get_simpleelem(elem))
               elem_is_blunt &= (n_elems_for_node[node] > 1) || (node < 4);

           elem_in_vacuum[elem] = elem_is_blunt ;
       }

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    bulk->copy_nodes(this->mesh);

    // Reserve memory for elements
    const int n_vacuum_elems = vector_sum(elem_in_vacuum);
    const int n_bulk_elems = n_elems - n_vacuum_elems;
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_bulk_elems);

    // Separate vacuum elements
    for (int i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i])
            vacuum->add_elem(mesh->get_simpleelem(i));
        else
            bulk->add_elem(mesh->get_simpleelem(i));

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

// Separate vacuum mesh from the union mesh by the element markers
// Handle edges that belong only to one element
const void Mesher::separate_meshes_vol3(Mesh* vacuum, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const int n_nodes = mesh->get_n_nodes();
    vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);

    // Get the number of elements the node is connected with
    vector<int> n_elems_for_node(n_nodes, 0);

    for (int elem = 0; elem < n_elems; ++elem)
        if (elem_in_vacuum[elem]) {
            for (int node : mesh->get_simpleelem(elem))
                n_elems_for_node[node]++;
        }

    // Assign elements that are sharp in vacuum to bulk where they are blunt
   for (int elem = 0; elem < n_elems; ++elem)
       if (elem_in_vacuum[elem]) {
           bool elem_is_blunt = true;

           // element is sharp if at least one of its nodes belongs to only one element
           // nodes 0-3 are the corners of surface, so skip them because their sharpness is OK
           for (int node : mesh->get_simpleelem(elem))
               elem_is_blunt &= (n_elems_for_node[node] > 1) || (node < 4);

           elem_in_vacuum[elem] = elem_is_blunt ;
       }

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    // Reserve memory for elements
    vacuum->init_elems(vector_sum(elem_in_vacuum));

    // Separate vacuum elements
    for (int i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i])
            vacuum->add_elem(mesh->get_simpleelem(i));

    vacuum->recalc(cmd);
}

/* Function to mark nodes with ray-triangle intersection technique
 *
 * When the ray from the node in z-direction crosses the surface, the node is located in material,
 * otherwise it's in vacuum.
 * The technique works perfectly with completely convex surface. The concaves give some false
 * nodes but because of their low spatial density they can be eliminated in post-processor.
 */
const void Mesher::mark_mesh(bool postprocess, double mean_thickness) {
    int node, elem, coord;
    const Vec3 ray_direction(0, 0, 1);
    Vec3 ray_origin;
    RaySurfaceIntersect rsi(mesh);
    rsi.generate_surf_faces();
    rsi.precompute_triangles(ray_direction);

    rsi.write_faces("output/faces_raysurf.vtk");

    // Calculate the min and max-s of nodes
    mesh->calc_statistics(0);
    // Mark faces by their location in relation to simulation cell edges
    mark_faces();

    // Reserve memory for the node markers
    mesh->init_nodemarkers(mesh->get_n_nodes());

    // Mark surface nodes by their known position in array
    for (node = mesh->indxs.surf_start; node <= mesh->indxs.surf_end; ++node)
        mesh->add_nodemarker(TYPES.SURFACE);

    // Mark bulk nodes by their known position in array
    for (node = mesh->indxs.bulk_start; node <= mesh->indxs.bulk_end; ++node)
        mesh->add_nodemarker(TYPES.BULK);

    // Mark vacuum nodes by their known position in array
    for (node = mesh->indxs.vacuum_start; node <= mesh->indxs.vacuum_end; ++node)
        mesh->add_nodemarker(TYPES.VACUUM);

    // Loop through all the nodes made by tetgen
    // and mark them by vertical-ray-intersects-surface-faces-technique
    for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node) {
        // Compose the ray_origin vector
        ray_origin = mesh->get_vec(node);

        // If ray intersects at least one surface triangle, the node
        // is considered to be located in bulk, otherwise it's located to vacuum
        if (rsi.ray_intersects_surface_fast(ray_origin, ray_direction))
            mesh->add_nodemarker(TYPES.BULK);
        else
            mesh->add_nodemarker(TYPES.VACUUM);
    }

    // Mark the elements by the node markers
    mark_elems();

    // Post process the nodes in shadow areas
    if (postprocess) post_process_marking();
}

// Force the bulk nodes in vacuum elements to become vacuum nodes
const void Mesher::post_process_marking() {
    const int n_elems = mesh->get_n_elems();

    bool node_changed = true;

    // Make as many re-check loops as necessary, but no more than 100
    for (int safety_cntr = 0; (safety_cntr < 100) && node_changed; ++safety_cntr) {
        // Post-process the wrongly marked nodes in vacuum
        node_changed = false;
        for (int elem = 0; elem < n_elems; ++elem) {
            if (mesh->get_elemmarker(elem) != TYPES.VACUUM) continue;

            // Force all the nodes in vacuum elements to be non-bulk ones
            for (int node : mesh->get_simpleelem(elem))
                if (mesh->get_nodemarker(node) == TYPES.BULK) {
                    mesh->set_nodemarker(node, TYPES.VACUUM);
                    node_changed = true;
                }
        }

        // If some of the nodes were changed,
        // mark all the surface and bulk elements again
        if (node_changed)
            for (int elem = 0; elem < n_elems; ++elem) {
                if (mesh->get_elemmarker(elem) == TYPES.VACUUM) continue;

                int location = 0;
                // Loop through all the nodes in element
                for (int node : mesh->get_simpleelem(elem)) {
                    int nodemarker = mesh->get_nodemarker(node);
                    // Encode element location into integer
                    if (nodemarker == TYPES.VACUUM)
                        location++;
                    else if (nodemarker == TYPES.BULK)
                        location--;
                }

                if (location > 0)
                    mesh->set_elemmarker(elem, TYPES.VACUUM);
                else if (location < 0)
                    mesh->set_elemmarker(elem, TYPES.BULK);
                else
                    mesh->set_elemmarker(elem, TYPES.SURFACE);
           }
    }

    // Among the other things calculate the number of atoms with given types
    mesh->calc_statistics(0);
    require(mesh->stat.n_bulk >= 4, "Nodemarker post-processor deleted the bulk atoms.\n"
            "Consider altering the surface refinement factor or disabling the post-processing.");
}

/* Function to mark nodes with ray-triangle intersection technique
 *
 * When the ray from the node in "random" direction crosses the surface 2*n times, the node
 * is located in vacuum, when it crosses the surface for 2*n+1, times, it's in material.
 * The technique works well with smooth surface. The non-smooth one gives for some ray direction
 * too many crossings. The workaround is to change the ray direction for the problematic nodes
 * to make sure the ray crosses only "quality" triangles.
 */
const void Mesher::mark_mesh_long() {
    int node, coord;
    Vec3 ray_origin, ray_direction(0.0, 0.0, 1.0);
    vector<int> markers(mesh->get_n_nodes(), TYPES.NONE);
    RaySurfaceIntersect rsi(mesh);
    rsi.generate_surf_faces();

    // Set of x & y directions for probing ray
    const double tilt_x[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!
    const double tilt_y[3] = { 0, 1.0, -1.0 }; // 0 must be the first entry!

    mesh->init_nodemarkers(mesh->get_n_nodes());

    // Mark surface nodes by their known position in array
    for (node = mesh->indxs.surf_start; node <= mesh->indxs.surf_end; ++node)
        mesh->add_nodemarker(TYPES.SURFACE);

    // Mark bulk nodes by their known position in array
    for (node = mesh->indxs.bulk_start; node <= mesh->indxs.bulk_end; ++node)
        mesh->add_nodemarker(TYPES.BULK);

    // Mark vacuum nodes by their known position in array
    for (node = mesh->indxs.vacuum_start; node <= mesh->indxs.vacuum_end; ++node)
        mesh->add_nodemarker(TYPES.VACUUM);

    bool all_nodes_done = false;

    // Try different ray directions until all the nodes have <= 1 intersection with surface faces
    for (const int tx : tilt_x) {
        if (all_nodes_done) break;

        ray_direction[0] = tx;
        for (const int ty : tilt_x) {
            if (all_nodes_done) break;

            ray_direction[1] = ty;
            ray_direction.normalize();

            rsi.precompute_triangles(ray_direction);

            all_nodes_done = true;

            // Loop through all the nodes inserted by tetgen
            for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node) {
                // If node is already checked, skip it
                if (markers[node] != TYPES.NONE) continue;

                // Compose the ray_origin vector
                ray_origin = mesh->get_vec(node);
                // Find how many times the ray from node crosses the surface
                int crossings = rsi.ray_intersects_surface_fast(ray_origin, ray_direction);

                if (crossings == 0)
                    markers[node] = TYPES.VACUUM;
                else if (crossings == 1)
                    markers[node] = TYPES.BULK;
                else
                    all_nodes_done = false;
            }
        }
    }

    // Insert found markers to mesh
    for (node = mesh->indxs.tetgen_start; node < mesh->get_n_nodes(); ++node)
        mesh->add_nodemarker(markers[node]);
}

// Mark elements by the node markers
const void Mesher::mark_elems() {
    const int n_elems = mesh->get_n_elems();

    mesh->init_elemmarkers(n_elems);

    // Loop through all the elements
    for (int elem = 0; elem < n_elems; ++elem) {
        int location = 0;

        // Loop through all the nodes in element
        for (int node : mesh->get_simpleelem(elem)) {
            int nodemarker = mesh->get_nodemarker(node);

            // Encode element location into integer
            if (nodemarker == TYPES.VACUUM)
                location++;
            else if (nodemarker == TYPES.BULK)
                location--;
        }

        // Element in vacuum is supposed not to have nodes in bulk
        if (location > 0)
            mesh->add_elemmarker(TYPES.VACUUM);
        // Element in bulk is supposed not to have nodes in vacuum
        else if (location < 0)
            mesh->add_elemmarker(TYPES.BULK);
        // Element in surface is supposed consist only of nodes in surface
        else
            mesh->add_elemmarker(TYPES.SURFACE);
    }
}

// Mark the boundary faces of mesh
const void Mesher::mark_faces() {
    const int eps = 0.1;
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);

    for (int i = 0; i < N; ++i) {
        Point3 centre = mesh->get_face_centre(i);

        if (on_boundary(centre.x, mesh->stat.xmin, eps))
            mesh->add_facemarker(TYPES.XMIN);
        else if (on_boundary(centre.x, mesh->stat.xmax, eps))
            mesh->add_facemarker(TYPES.XMAX);
        else if (on_boundary(centre.y, mesh->stat.ymin, eps))
            mesh->add_facemarker(TYPES.YMIN);
        else if (on_boundary(centre.y, mesh->stat.ymax, eps))
            mesh->add_facemarker(TYPES.YMAX);
        else if (on_boundary(centre.z, mesh->stat.zmin, eps))
            mesh->add_facemarker(TYPES.ZMIN);
        else if (on_boundary(centre.z, mesh->stat.zmax, eps))
            mesh->add_facemarker(TYPES.ZMAX);
        else
            mesh->add_facemarker(TYPES.SURFACE);
    }
}

} /* namespace femocs */
