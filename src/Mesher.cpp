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
//    faces.reserve(0);
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
    const int n_faces = mesh->get_n_faces();// faces.size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        SimpleFace sface = mesh->get_simpleface(face);
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
    const int n_faces = mesh->get_n_faces();
    int crossings = 0;

    vector<bool> node_passed(mesh->get_n_nodes(), false);

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        if (is_parallel[face]) continue;

        SimpleFace sface = mesh->get_simpleface(face);

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
const bool RaySurfaceIntersect::ray_intersects_surface_fast(const Vec3 &origin, const Vec3 &direction) {
    const int n_faces = mesh->get_n_faces();

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face) {
        if (ray_intersects_triangle(origin, direction, face)) return true;
    }

    return false;
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

    // Memorize the positions of different kind of nodes to speed up later calculations
    mesh->indxs.surf_start = 0;
    mesh->indxs.surf_end = mesh->indxs.surf_start + n_surf - 1;
    mesh->indxs.bulk_start = mesh->indxs.surf_end + 1;
    mesh->indxs.bulk_end = mesh->indxs.bulk_start + n_bulk - 1;
    mesh->indxs.vacuum_start = mesh->indxs.bulk_end + 1;
    mesh->indxs.vacuum_end = mesh->indxs.vacuum_start + n_vacuum - 1;
    mesh->indxs.tetgen_start = mesh->indxs.vacuum_end + 1;

    mesh->recalc("Q", cmd);
}

// Generate manually edges and surface faces
const void Mesher::generate_mesh_appendices() {
    // Generate edges from elements
//    generate_edges();
    // Generate surface faces from elements
    generate_surf_faces();
}

// Function to generate edges from the elements
// Overlapping edges are not cleaned to make it easier to match them with element
const void Mesher::generate_edges() {
    const int n_elems = mesh->get_n_elems();
    const int n_edges = mesh->n_edges_per_elem * n_elems;

    mesh->init_edges(n_edges);

    // Generate edges from all the elements and leave
    for (int elem = 0; elem < n_elems; ++elem) {
        SimpleElement selem = mesh->get_simpleelem(elem);
        cout << selem.edge(0);
        for (unsigned int ed = 0; ed < mesh->n_edges_per_elem; ++ed)
            mesh->add_edge(selem.edge(ed));
    }
}

// Function to manually generate surface faces from elements and surface nodes
// Overlapping faces are not cleaned for computational efficiency purposes
const void Mesher::generate_surf_faces() {
    const int n_elems = mesh->get_n_elems();
    const int max_surf_indx = mesh->indxs.surf_end;
    const int* elems = mesh->get_elems();

    // booleans showing whether element i has exactly one face on the surface or not
    vector<bool> elem_in_surface(n_elems);
    // booleans showing whether node i is on the surface or not
    vector<bool> surf_locs;

    // Mark the elements that have exactly one face on the surface
    for (int elem = 0; elem < n_elems; ++elem) {
        surf_locs = mesh->get_simpleelem(elem) <= max_surf_indx;
        elem_in_surface[elem] = (vector_sum(surf_locs) == 3);
    }
    // Reserve memory for surface faces
    mesh->init_faces( 2 + vector_sum(elem_in_surface) );

    // Make two big faces that pass the xy-min-max corners of surface
    mesh->add_face(SimpleFace(0, 1, 2));
    mesh->add_face(SimpleFace(0, 2, 3));

    // Generate the faces that separate material and vacuum
    // The faces are taken from the elements that have exactly one face on the surface
    for (int elem = 0; elem < n_elems; ++elem)
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
            mesh->add_face(SimpleFace(n0, n1, n2));
        }
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
            vacuum->add_elem(vacm_map[se[0]], vacm_map[se[1]], vacm_map[se[2]], vacm_map[se[3]]);
        else
            bulk->add_elem(bulk_map[se[0]], bulk_map[se[1]], bulk_map[se[2]], bulk_map[se[3]]);
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
            vacuum->add_elem(vacm_map[se[0]], vacm_map[se[1]], vacm_map[se[2]], vacm_map[se[3]]);
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

// Separate vacuum and bulk mesh from the union mesh by the element markers
// Handle edges that belong only to one element
const void Mesher::separate_meshes_vol2(Mesh* vacuum, Mesh* bulk, const string& cmd) {
    const int n_elems = mesh->get_n_elems();
    const int n_nodes = mesh->get_n_nodes();
//    vector<bool> elem_in_vacuum = vector_not(mesh->get_elemmarkers(), TYPES.BULK);
    vector<bool> elem_in_vacuum = vector_equal(mesh->get_elemmarkers(), TYPES.VACUUM);

    // Mark elements that are on the perimeter of simulation box
//    vector<bool> elem_on_perim(n_elems);
//    for (int i = 0; i < n_elems; ++i) {
//        int I = i * mesh->n_edges_per_elem;
//        bool on_edge = false;
//        for (int j = 0; j < mesh->n_edges_per_elem; ++j)
//            on_edge |= mesh->get_edgemarker(I+j) == TYPES.PERIMETER;
//        elem_on_perim[i] = on_edge;
//    }
//
//    swap_sharp_elements(elem_in_vacuum, elem_on_perim);
//    swap_sharp_elements(elem_in_vacuum, elem_on_perim);
//    swap_sharp_elements(elem_in_vacuum, elem_on_perim);

    // Copy the non-bulk nodes from input mesh without modification
    vacuum->copy_nodes(this->mesh);
    bulk->copy_nodes(this->mesh);

    // Reserve memory for elements
    const int n_elems_vacuum = vector_sum(elem_in_vacuum);
    const int n_elems_bulk = n_elems - n_elems_vacuum;
    vacuum->init_elems(n_elems_vacuum);
    bulk->init_elems(n_elems_bulk);

    // Separate vacuum and bulk elements
    for (int elem = 0; elem < n_elems; ++elem)
        if (elem_in_vacuum[elem])
            vacuum->add_elem(mesh->get_simpleelem(elem));
        else
            bulk->add_elem(mesh->get_simpleelem(elem));

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

const void Mesher::swap_sharp_elements(vector<bool> &elem_in_vacuum, vector<bool> &elem_on_perim) {
    const int n_elems = mesh->get_n_elems();

    // Find the number of neighbours for vacuum elements
    vector<int> nn_vacuum(n_elems, 0);
    for (int elem = 0; elem < n_elems; ++elem) {
        if (elem_in_vacuum[elem])
            for (int node : mesh->get_elem_neighbours(elem))
                if (node >= 0 && elem_in_vacuum[node])
                    nn_vacuum[elem]++;
    }

//    // Find the number of neighbours for bulk elements
//    vector<int> nn_bulk(n_elems, 0);
//    for (int elem = 0; elem < n_elems; ++elem) {
//        if (!elem_in_vacuum[elem])
//            for (int node : mesh->get_elem_neighbours(elem))
//                if (node >= 0 && !elem_in_vacuum[node])
//                    nn_bulk[elem]++;
//    }

    for (int elem = 0; elem < n_elems; ++elem) {
        if (elem_on_perim[elem])
            continue;
        else if (nn_vacuum[elem] < 3)
            elem_in_vacuum[elem] = false;
//        else if (nn_bulk[elem] < 3)
//            elem_in_vacuum[elem] = true;
    }
}

// Mark mesh nodes, edges, faces and elements
const bool Mesher::mark_mesh(const bool postprocess) {
    bool marking_successful = true;

    // Mark nodes with ray-triangle intersection technique
    mark_nodes();

    // Mark the elements by the node markers
    mark_elems();

    // Post process the nodes in shadow areas
    if (postprocess) marking_successful = post_process_marking();

    if (!marking_successful)
        return false;

//    // Mark nodes on simulation cell perimeter
//    remark_perimeter_nodes();
//
//    // Mark edges by the node markers
//    mark_edges();

    // Mark faces on simulation cell edges
//    mark_faces;

    return marking_successful;
}

// Force the bulk nodes in vacuum elements to become vacuum nodes
const bool Mesher::post_process_marking() {
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

        // If some of the nodes were changed, mark all the surface and bulk elements again
        if (node_changed) remark_elems(TYPES.VACUUM);
    }

    // Among the other things calculate the number of atoms with given types
    mesh->calc_statistics(0);
    expect(mesh->stat.n_bulk > 4, "Nodemarker post-processor deleted the bulk atoms.\n"
            "Consider altering the surface refinement factor or disabling the post-processing.");

    // Return value indicates whether the new system is valid or not
    return mesh->stat.n_bulk > 4;
}

/* Function to mark nodes with ray-triangle intersection technique
 *
 * When the ray from the node in z-direction crosses the surface, the node is located in material,
 * otherwise it's in vacuum.
 * The technique works perfectly with completely convex surface. The concaves give some false
 * nodes but because of their low spatial density they can be eliminated in post-processor.
 */
const void Mesher::mark_nodes() {
    int node, elem, coord;
    const Vec3 ray_direction(0, 0, 1);
    RaySurfaceIntersect rsi(mesh);
    rsi.precompute_triangles(ray_direction);

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
        Vec3 ray_origin = mesh->get_vec(node);

        // If ray intersects at least one surface triangle, the node
        // is considered to be located in bulk, otherwise it's located in vacuum
        if (rsi.ray_intersects_surface_fast(ray_origin, ray_direction))
            mesh->add_nodemarker(TYPES.BULK);
        else
            mesh->add_nodemarker(TYPES.VACUUM);
    }
}

// Remark vacuum nodes on simulation cell perimeter
const void Mesher::remark_perimeter_nodes() {
    const int eps = 0.1;
    const int n_nodes = mesh->get_n_nodes();

    for (int node = 0; node < n_nodes; ++node) {
        Point3 point = mesh->get_node(node);
        // Point on xmin or xmax of simubox?
        const bool bound_x = on_boundary(point.x, mesh->stat.xmin, mesh->stat.xmax, eps);
        // Point on ymin or ymax of simubox?
        const bool bound_y = on_boundary(point.y, mesh->stat.ymin, mesh->stat.ymax, eps);
        // Point on zmax of simubox?
        const bool bound_z = on_boundary(point.z, mesh->stat.zmax, eps);
        // Point on the perimeter of surface?
        const bool bound_surf = (mesh->get_nodemarker(node) == TYPES.SURFACE) && (bound_x || bound_y);
        // Point not in the bulk area?
        const bool not_bulk = mesh->get_nodemarker(node) != TYPES.BULK;

        if ( not_bulk && ((bound_x && bound_y) || (bound_x && bound_z) || (bound_y && bound_z) || bound_surf) )
            mesh->set_nodemarker(node, TYPES.PERIMETER);
    }
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

// Remark elements skipping the ones marked as skip_type
const void Mesher::remark_elems(const int skip_type) {
    // Loop through all the non-vacuum elements
    for (int elem = 0; elem < mesh->get_n_elems(); ++elem) {
        if (mesh->get_elemmarker(elem) == skip_type) continue;

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

// Mark the edges on the simulation cell perimeter by the node markers
const void Mesher::mark_edges() {
    const int n_edges = mesh->get_n_edges();
    mesh->init_edgemarkers(n_edges);

    for (int edge = 0; edge < n_edges; ++edge) {
        SimpleEdge sedge = mesh->get_simpleedge(edge);
        const int m1 = mesh->get_nodemarker(sedge[0]);
        const int m2 = mesh->get_nodemarker(sedge[1]);

        if (m1 == TYPES.PERIMETER && m2 == TYPES.PERIMETER)
            mesh->add_edgemarker(TYPES.PERIMETER);
        else
            mesh->add_edgemarker(TYPES.NONE);
    }
}

// Mark the boundary faces of mesh
const void Mesher::mark_faces() {
    const int eps = 0.1;
    const int n_faces = mesh->get_n_faces();

    mesh->init_facemarkers(n_faces);

    for (int i = 0; i < n_faces; ++i) {
        Point3 centre = mesh->get_face_centroid(i);

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
