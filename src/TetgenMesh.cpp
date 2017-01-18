/*
 * TetgenMesh.cpp
 *
 *  Created on: 3.10.2016
 *      Author: veske
 */

#include "TetgenMesh.h"
#include "Tethex.h"
#include <fstream>

using namespace std;
namespace femocs {

RaySurfaceIntersect::RaySurfaceIntersect(TetgenMesh* mesh) : mesh(mesh) {
    reserve(0);
}

// Reserve memory for precompute data
const void RaySurfaceIntersect::reserve(const int n) {
    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear(); pvec.reserve(n);
    is_parallel.clear(); is_parallel.reserve(n);
}

// Precompute the data needed to execute the Moller-Trumbore algorithm
const void RaySurfaceIntersect::precompute_triangles(const Vec3 &direction) {
    const int n_faces = mesh->faces.size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Loop through all the faces
    for (SimpleFace sface : mesh->faces) {
        Vec3 v0 = mesh->nodes.get_vec(sface[0]);
        Vec3 v1 = mesh->nodes.get_vec(sface[1]);
        Vec3 v2 = mesh->nodes.get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = direction.crossProduct(e2);
        double det = e1.dotProduct(pv);
        double i_det = 1.0 / det;

        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        vert0.push_back(v0);
        pvec.push_back(pv * i_det);
        is_parallel.push_back(fabs(det) < epsilon);
    }
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

// Determine whether the ray intersects with any of the triangles on the surface mesh
const bool RaySurfaceIntersect::ray_intersects_surface(const Vec3 &origin, const Vec3 &direction) {
    const int n_faces = mesh->faces.size();

    // Loop through all the faces
    for (int face = 0; face < n_faces; ++face)
        if (ray_intersects_triangle(origin, direction, face)) return true;

    return false;
}

// Group hexahedra around tetrahedra nodes
const void TetgenMesh::group_hexs() {
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // find which hexahedra correspond to which tetrahedral node
    // hexahedra with the same tetrahedral node form the voronoi cell of that node
    for (int i = 0; i < hexahedra.size(); ++i) {
        for (int node : hexahedra[i])
            if (node >= node_min && node <= node_max) {
                hexahedra.set_marker(i, node);
                continue;
            }
    }
}

// Group hexahedra around tetrahedra nodes
const void TetgenMesh::get_voronoi_cells() {
    vector<vector<unsigned int>> voronoi_cells(nodes.stat.n_tetnode);
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // enumerate tetrahedral nodes
    for (int i = 0; i < nodes.size(); ++i)
        if (nodes.get_marker(i) == 1)
            nodes.set_marker(i, -i);

    // find the voronoi cell nodes for the tetrahedral nodes
    for (int i = 0; i < hexahedra.size(); ++i) {
        for (unsigned int node : hexahedra[i])
            if (nodes.get_marker(node) == 4) {
                const int tetnode = hexahedra.get_marker(i);
                expect(tetnode >= node_min &&  tetnode <= node_max, "Illegal node detected: " + to_string(tetnode));
                voronoi_cells[tetnode].push_back(node);
//                nodes.set_marker(node, tetnode);
                continue;
            }
    }

    // mark the voronoi cell nodes
    for (int i = 0; i < voronoi_cells.size(); ++i) {
        cout << i << " \t" << voronoi_cells[i].size() << ": ";
        for (int v : voronoi_cells[i]) {
            nodes.set_marker(v, i);
            cout << v << " ";
        }
        cout << endl;
    }
}

// Function to generate simple mesh that consists of one tetrahedron
const bool TetgenMesh::generate_simple() {
    const int n_nodes = elems.DIM;
    const int n_edges = n_edges_per_elem;
    const int n_faces = n_faces_per_elem;
    const int n_elems = 1;

    nodes.init(n_nodes);
    nodes.append(Point3(1.0, 0.0, 0.7));
    nodes.append(Point3(-1.0, 0.0, 0.7));
    nodes.append(Point3(0.0, 1.0, -0.7));
    nodes.append(Point3(0.0, -1.0, -0.7));

    faces.init(n_faces);
    faces.append(SimpleFace(0, 1, 3));
    faces.append(SimpleFace(1, 2, 3));
    faces.append(SimpleFace(2, 0, 3));
    faces.append(SimpleFace(0, 1, 2));

    edges.init(n_edges);
    edges.append(SimpleEdge(0, 1));
    edges.append(SimpleEdge(0, 2));
    edges.append(SimpleEdge(0, 3));
    edges.append(SimpleEdge(1, 2));
    edges.append(SimpleEdge(1, 3));
    edges.append(SimpleEdge(2, 3));

    elems.init(n_elems);
    elems.append(SimpleElement(0, 1, 2, 3));

    return recalc("rQ");
}

// Copy mesh from input to output without modification
const bool TetgenMesh::recalc() {
    nodes.recalc();
    edges.recalc();
    faces.recalc();
    elems.recalc();
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
const bool TetgenMesh::recalc(const string& cmd) {
    try {
        tetrahedralize(const_cast<char*>(cmd.c_str()), &tetIOin, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        faces.set_counter(tetIOout.numberoftrifaces);
        elems.set_counter(tetIOout.numberoftetrahedra);
    }
    catch (int e) { return 1; }
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
const bool TetgenMesh::recalc(const string& cmd1, const string& cmd2) {
    try {
        tetgenio tetIOtemp;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp);
        tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        faces.set_counter(tetIOout.numberoftrifaces);
        elems.set_counter(tetIOout.numberoftetrahedra);
    }
    catch (int e) { return 1; }
    return 0;
}

// Write mesh into files with Tetgen functions
const bool TetgenMesh::write_tetgen(const string& file_name) {
    const string cmd = "Q";
    tetgenbehavior tetgenbeh;

    try {
        tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
        for (int i = 0; i < file_name.size(); ++i)
            tetgenbeh.outfilename[i] = file_name[i];

        tetrahedralize(&tetgenbeh, &tetIOout, NULL);
    } catch (int e) { return 1; }

    return 0;
}

// Function to generate mesh from surface, bulk and vacuum atoms
const bool TetgenMesh::generate(const Media& bulk, const Media& surf, const Media& vacuum, const string& cmd) {
    const int n_bulk = bulk.get_n_atoms();
    const int n_surf = surf.get_n_atoms();
    const int n_vacuum = vacuum.get_n_atoms();

    nodes.init(n_bulk + n_surf + n_vacuum);

    // Add surface atoms first,...
    for (int i = 0; i < n_surf; ++i)
        nodes.append(surf.get_point(i));

    // ... bulk atoms second,...
    for (int i = 0; i < n_bulk; ++i)
        nodes.append(bulk.get_point(i));

    // ... and vacuum atoms last
    for (int i = 0; i < n_vacuum; ++i)
        nodes.append(vacuum.get_point(i));

    // Memorize the positions of different types of nodes to speed up later calculations
    nodes.save_indices(n_surf, n_bulk, n_vacuum);

    return recalc("Q", cmd);
}

const bool TetgenMesh::generate_hexs() {
    tethex::Mesh hexmesh;
    hexmesh.read_femocs(this);
    hexmesh.convert();
    hexmesh.export_femocs(this);

    return 0;
}

// Generate manually edges and surface faces
const bool TetgenMesh::generate_appendices() {
    // Generate edges from elements
//    generate_edges();
    // Generate surface faces from elements
    generate_surf_faces();
    return 0;
}

// Function to generate edges from the elements
// Overlapping edges are not cleaned to make it easier to match them with element
const void TetgenMesh::generate_edges() {
    const int n_elems = elems.size();
    const int n_edges = n_edges_per_elem * n_elems;

    edges.init(n_edges);

    // Generate edges from all the elements and leave
    for (SimpleElement selem : elems) {
        for (int e = 0; e < n_edges_per_elem; ++e)
            edges.append(selem.edge(e));
    }
}

// Function to manually generate surface faces from elements and surface nodes
// Overlapping faces are not cleaned for computational efficiency purposes
const void TetgenMesh::generate_surf_faces() {
    const int n_elems = elems.size();
    const int max_surf_indx = nodes.indxs.surf_end;

    // booleans showing whether element i has exactly one face on the surface or not
    vector<bool> elem_on_surface; elem_on_surface.reserve(n_elems);
    // booleans showing whether node i is on the surface or not
    vector<bool> surf_locs;

    // Mark the elements that have exactly one face on the surface
    for (SimpleElement selem : elems) {
        surf_locs = (selem <= max_surf_indx);
        elem_on_surface.push_back(vector_sum(surf_locs) == 3);
    }

    // Reserve memory for surface faces
    faces.init( vector_sum(elem_on_surface) + 2 );

    // Make two big faces that pass the xy-min-max corners of surface
    faces.append(SimpleFace(0, 1, 2));
    faces.append(SimpleFace(0, 2, 3));

    // Generate the faces that separate material and vacuum
    // The faces are taken from the elements that have exactly one face on the surface
    for (int el = 0; el < n_elems; ++el)
        if (elem_on_surface[el]) {
            SimpleElement elem = elems[el];

            // Find the indices of nodes that are on the surface
            surf_locs = (elem <= max_surf_indx);

            /* The possible combinations of surf_locs and n0,n1,n2:
             * surf_locs: 1110   1101   1011   0111
             *        n0: elem0  elem0  elem0  elem1
             *        n1: elem1  elem1  elem2  elem2
             *        n2: elem2  elem3  elem3  elem3   */
            int n0 = surf_locs[0] * elem[0] + (!surf_locs[0]) * elem[1];
            int n1 = (surf_locs[0] & surf_locs[1]) * elem[1] + (surf_locs[2] & surf_locs[3]) * elem[2];
            int n2 = (!surf_locs[3]) * elem[2] + surf_locs[3] * elem[3];
            faces.append(SimpleFace(n0, n1, n2));
        }
}

// Separate vacuum and bulk mesh from the union mesh by the element markers
const bool TetgenMesh::separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd) {
    vector<bool> tet_mask = vector_equal(elems.get_markers(), TYPES.VACUUM);
    vector<bool> hex_mask = vector_greater_equal(hexahedra.get_markers(), 0);

    // Transfer vacuum nodes, tetrahedra, hexahedra and their markers
    vacuum.nodes.copy(this->nodes);
    vacuum.nodes.copy_markers(this->nodes);
    vacuum.elems.copy(this->elems, tet_mask);
    vacuum.elems.copy_markers(this->elems, tet_mask);
    vacuum.hexahedra.copy(this->hexahedra, hex_mask);
    vacuum.hexahedra.copy_markers(this->hexahedra, hex_mask);

    tet_mask.flip();
    hex_mask.flip();

    // Transfer bulk nodes, tetrahedra, hexahedra and their markers
    bulk.nodes.copy(this->nodes);
    bulk.nodes.copy_markers(this->nodes);
    bulk.elems.copy(this->elems, tet_mask);
    bulk.elems.copy_markers(this->elems, tet_mask);
    bulk.hexahedra.copy(this->hexahedra, hex_mask);
    bulk.hexahedra.copy_markers(this->hexahedra, hex_mask);

    return vacuum.recalc(cmd) ||  bulk.recalc(cmd);
}

const bool TetgenMesh::smoothen(double radius, double smooth_factor, double r_cut) {
    const int n_atoms = nodes.size();

    // Find atoms in vacuum-bulk boundary
    vector<bool> hot_node = vector_equal(nodes.get_markers(), TYPES.SURFACE);

    // Turn atoms on boundary into surface
    Media surf(vector_sum(hot_node));

    int cntr = 0;
    for (int i = 0; i < n_atoms; ++i)
        if (hot_node[i])
            surf.add_atom( nodes[i] );

    // Smoothen the surface
    surf.smoothen(radius, smooth_factor, r_cut);

    // Write smoothed vertices back to mesh
    int j = 0;
    for(int i = 0; i < n_atoms; ++i)
        if(hot_node[i])
            nodes.set_node(i, surf.get_point(j++));

    return 0;
}

// Mark mesh nodes, edges, faces and elements
const bool TetgenMesh::mark_mesh(const bool postprocess) {
//    mark_elems_byrsi();
//    mark_nodes();
//    return true;

    bool fail = 0;

    // Mark nodes with ray-triangle intersection technique
    mark_nodes();

    // Mark the elements by the node markers
    mark_elems_vol2();

    // Post process the nodes in shadow areas
    if (postprocess)
        fail = post_process_marking();

    if (fail) return 1;

//    // Mark nodes on simulation cell perimeter
//    remark_perimeter_nodes();
//
//    // Mark edges by the node markers
//    mark_edges();

    // Mark faces on simulation cell edges
//    mark_faces;

    return 0;
}

// Force the bulk nodes in vacuum elements to become vacuum nodes
const bool TetgenMesh::post_process_marking() {
    const int n_elems = elems.size();

    bool node_changed = true;

    // Post-process the wrongly marked nodes in vacuum
    while (node_changed) {
        node_changed = false;
        for (int elem = 0; elem < n_elems; ++elem) {
            if (elems.get_marker(elem) != TYPES.VACUUM) continue;

            // Force all the nodes in vacuum elements to be non-bulk ones
            for (int node : elems[elem])
                if (nodes.get_marker(node) == TYPES.BULK) {
                    nodes.set_marker(node, TYPES.VACUUM);
                    node_changed = true;
                }
        }

        // If some of the nodes were changed, mark all the surface and bulk elements again
        if (node_changed) remark_elems(TYPES.VACUUM);
    }

    // Among the other things calculate the number of atoms with given types
    nodes.calc_statistics();
    expect(nodes.stat.n_bulk > 4, "Nodemarker post-processor deleted the bulk atoms.\n"
            "Consider altering the surface refinement factor or disabling the post-processing.");

    // Return value indicates whether the new system is invalid or not
    return nodes.stat.n_bulk <= 4;
}

/* Function to mark nodes with ray-triangle intersection technique
 *
 * When the ray from the node in z-direction crosses the surface, the node is located in material,
 * otherwise it's in vacuum.
 * The technique works perfectly with completely convex surface. The concaves give some false
 * nodes but because of their low spatial density they can be eliminated in post-processor.
 */
const void TetgenMesh::mark_nodes() {
    int node, elem, coord;
    const Vec3 ray_direction(0, 0, 1);
    RaySurfaceIntersect rsi(this);
    rsi.precompute_triangles(ray_direction);

    // Calculate the min and max-s of nodes
    nodes.calc_statistics();
    // Mark faces by their location in relation to simulation cell edges
    mark_faces();

    // Reserve memory for the node markers
    nodes.init_markers(nodes.size());

    // Mark surface nodes by their known position in array
    for (node = nodes.indxs.surf_start; node <= nodes.indxs.surf_end; ++node)
        nodes.append_marker(TYPES.SURFACE);

    // Mark bulk nodes by their known position in array
    for (node = nodes.indxs.bulk_start; node <= nodes.indxs.bulk_end; ++node)
        nodes.append_marker(TYPES.BULK);

    // Mark vacuum nodes by their known position in array
    for (node = nodes.indxs.vacuum_start; node <= nodes.indxs.vacuum_end; ++node)
        nodes.append_marker(TYPES.VACUUM);

    // Loop through all the nodes made by tetgen
    // and mark them by vertical-ray-intersects-surface-faces-technique
    for (node = nodes.indxs.tetgen_start; node < nodes.size(); ++node) {
        // Compose the ray_origin vector
        Vec3 ray_origin = nodes.get_vec(node);

        // If ray intersects at least one surface triangle, the node
        // is considered to be located in bulk, otherwise it's located in vacuum
        if (rsi.ray_intersects_surface(ray_origin, ray_direction))
            nodes.append_marker(TYPES.BULK);
        else
            nodes.append_marker(TYPES.VACUUM);
    }
}

const void TetgenMesh::mark_elems_byrsi() {
    int node, elem, coord;
    const Vec3 ray_direction(0, 0, 1);
    RaySurfaceIntersect rsi(this);
    rsi.precompute_triangles(ray_direction);

    // Reserve memory for the element markers
    elems.init_markers(elems.size());

    // Loop through all the nodes made by tetgen
    // and mark them by vertical-ray-intersects-surface-faces-technique
    for (elem = 0; elem < elems.size(); ++elem) {
        // if all the nodes of element are located on surface,
        // the element is considered to be located on surface
        SimpleElement selem = elems[elem];
        if ( vector_sum(selem < nodes.indxs.surf_end) == selem.size() ) {
            elems.append_marker(TYPES.SURFACE);
            continue;
        }

        // Compose the ray_origin vector
        Vec3 ray_origin(elems.get_centroid(elem));

        // If ray intersects at least one surface triangle, the element
        // is considered to be located in bulk, otherwise it's located in vacuum
        if (rsi.ray_intersects_surface(ray_origin, ray_direction))
            elems.append_marker(TYPES.BULK);
        else
            elems.append_marker(TYPES.VACUUM);
    }
}

// Remark vacuum nodes on simulation cell perimeter
const void TetgenMesh::remark_perimeter_nodes() {
    const int eps = 0.1;
    const int n_nodes = nodes.size();

    for (int node = 0; node < n_nodes; ++node) {
        Point3 point = nodes[node];
        // Point on xmin or xmax of simubox?
        const bool bound_x = on_boundary(point.x, nodes.stat.xmin, nodes.stat.xmax, eps);
        // Point on ymin or ymax of simubox?
        const bool bound_y = on_boundary(point.y, nodes.stat.ymin, nodes.stat.ymax, eps);
        // Point on zmax of simubox?
        const bool bound_z = on_boundary(point.z, nodes.stat.zmax, eps);
        // Point on the perimeter of surface?
        const bool bound_surf = (nodes.get_marker(node) == TYPES.SURFACE) && (bound_x || bound_y);
        // Point not in the bulk area?
        const bool not_bulk = nodes.get_marker(node) != TYPES.BULK;

        if ( not_bulk && ((bound_x && bound_y) || (bound_x && bound_z) || (bound_y && bound_z) || bound_surf) )
            nodes.set_marker(node, TYPES.PERIMETER);
    }
}

// Mark elements by the node markers
const void TetgenMesh::mark_elems() {
    const int n_elems = elems.size();

    elems.init_markers(n_elems);

    // Loop through all the elements
    for (SimpleElement elem : elems) {
        int location = 0;

        // Loop through all the nodes in element
        for (int node : elem) {
            const int nodemarker = nodes.get_marker(node);

            // Encode element location into integer
            if (nodemarker == TYPES.VACUUM)
                location++;
            else if (nodemarker == TYPES.BULK)
                location--;
        }

        // Element in vacuum is supposed not to have nodes in bulk
        if (location > 0)
            elems.append_marker(TYPES.VACUUM);
        // Element in bulk is supposed not to have nodes in vacuum
        else if (location < 0)
            elems.append_marker(TYPES.BULK);
        // Element in surface is supposed consist only of nodes in surface
        else
            elems.append_marker(TYPES.SURFACE);
    }
}

const void TetgenMesh::mark_elems_vol2() {
    const int n_elems = elems.size();

    elems.init_markers(n_elems);

    // Loop through all the elements
    for (SimpleElement elem : elems) {
        int n_vacuum = 0;
        int n_bulk = 0;

        // Loop through all the nodes in element
        for (int node : elem) {
            const int nodemarker = nodes.get_marker(node);

            // Encode element location into integer
            if (nodemarker == TYPES.VACUUM)
                n_vacuum++;
            else if (nodemarker == TYPES.BULK)
                n_bulk++;
        }

        // Element in vacuum is supposed not to have nodes in bulk
        if (n_vacuum > 0 && n_bulk == 0)
            elems.append_marker(TYPES.VACUUM);
        // Element in bulk is supposed not to have nodes in vacuum
        else if (n_bulk > 0 && n_vacuum == 0)
            elems.append_marker(TYPES.BULK);
        // Element in surface is supposed consist only of nodes in surface
        else if (n_bulk == 0 && n_vacuum == 0)
            elems.append_marker(TYPES.SURFACE);
        else
            elems.append_marker(TYPES.NONE);
    }
}

// Remark elements skipping the ones marked as skip_type
const void TetgenMesh::remark_elems(const int skip_type) {
    // Loop through all the non-vacuum elements
    for (int elem = 0; elem < elems.size(); ++elem) {
        if (elems.get_marker(elem) == skip_type) continue;

        int location = 0;
        // Loop through all the nodes in element
        for (int node : elems[elem]) {
            int nodemarker = nodes.get_marker(node);
            // Encode element location into integer
            if (nodemarker == TYPES.VACUUM)
                location++;
            else if (nodemarker == TYPES.BULK)
                location--;
        }

        if (location > 0)
            elems.set_marker(elem, TYPES.VACUUM);
        else if (location < 0)
            elems.set_marker(elem, TYPES.BULK);
        else
            elems.set_marker(elem, TYPES.SURFACE);
     }
}

// Mark the edges on the simulation cell perimeter by the node markers
const void TetgenMesh::mark_edges() {
    const int n_edges = edges.size();
    edges.init_markers(n_edges);

    for (int edge = 0; edge < n_edges; ++edge) {
        SimpleEdge sedge = edges[edge];
        const int m1 = nodes.get_marker(sedge[0]);
        const int m2 = nodes.get_marker(sedge[1]);

        if (m1 == TYPES.PERIMETER && m2 == TYPES.PERIMETER)
            edges.append_marker(TYPES.PERIMETER);
        else
            edges.append_marker(TYPES.NONE);
    }
}

// Mark the boundary faces of mesh
const void TetgenMesh::mark_faces() {
    const int eps = 0.1;
    const int n_faces = faces.size();

    faces.init_markers(n_faces);
    nodes.calc_statistics();

    for (int i = 0; i < n_faces; ++i) {
        Point3 centre = faces.get_centroid(i);

        if (on_boundary(centre.x, nodes.stat.xmin, eps))
            faces.append_marker(TYPES.XMIN);
        else if (on_boundary(centre.x, nodes.stat.xmax, eps))
            faces.append_marker(TYPES.XMAX);
        else if (on_boundary(centre.y, nodes.stat.ymin, eps))
            faces.append_marker(TYPES.YMIN);
        else if (on_boundary(centre.y, nodes.stat.ymax, eps))
            faces.append_marker(TYPES.YMAX);
        else if (on_boundary(centre.z, nodes.stat.zmin, eps))
            faces.append_marker(TYPES.ZMIN);
        else if (on_boundary(centre.z, nodes.stat.zmax, eps))
            faces.append_marker(TYPES.ZMAX);
        else
            faces.append_marker(TYPES.SURFACE);
    }
}
} /* namespace femocs */
