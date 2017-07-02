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

RaySurfaceIntersect::RaySurfaceIntersect(const TetgenMesh* m) : mesh(m) {
    reserve(0);
}

// Reserve memory for precompute data
void RaySurfaceIntersect::reserve(const int n) {
    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear(); pvec.reserve(n);
    is_parallel.clear(); is_parallel.reserve(n);
}

// Precompute the data needed to calculate the distance of points from surface in given direction
void RaySurfaceIntersect::precompute_triangles(const Vec3 &direction) {
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

// Precompute the data needed to calculate the distance of points from surface in the direction of triangle norms
void RaySurfaceIntersect::precompute_triangles() {
    const int n_faces = mesh->faces.size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace sface = mesh->faces[i];
        Vec3 v0 = mesh->nodes.get_vec(sface[0]);
        Vec3 v1 = mesh->nodes.get_vec(sface[1]);
        Vec3 v2 = mesh->nodes.get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = mesh->faces.get_norm(i).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);
    }
}

// Moller-Trumbore algorithm to find whether the ray and the triangle intersect or not
bool RaySurfaceIntersect::ray_intersects_triangle(const Vec3 &origin, const Vec3 &direction,
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

// Moller-Trumbore algorithm to get the distance of point from the triangle in the direction of triangle norm
double RaySurfaceIntersect::distance_from_triangle(const Vec3 &point, const int face) {
    double u, v;

    Vec3 tvec = point - vert0[face];
    u = tvec.dotProduct(pvec[face]);
    if (u < zero || u > one) return 1e100;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    v = mesh->faces.get_norm(face).dotProduct(qvec);
    if (v < zero || u + v > one) return 1e100; // Check second & third barycentric coordinate

    // return the distance from point to triangle
//    return fabs(edge2[face].dotProduct(qvec));
    return edge2[face].dotProduct(qvec);
}

// Moller-Trumbore algorithm to get the distance of point from the triangle in given direction
double RaySurfaceIntersect::distance_from_triangle(const Vec3 &point, const Vec3 &direction, const int face) {
    const SimpleFace sface = mesh->faces[face];
    const Vec3 v0 = mesh->nodes.get_vec(sface[0]);
    const Vec3 v1 = mesh->nodes.get_vec(sface[1]);
    const Vec3 v2 = mesh->nodes.get_vec(sface[2]);

    double u, v, det, invDet;

    Vec3 e1 = v1 - v0;
    Vec3 e2 = v2 - v0;
    Vec3 pvec = direction.crossProduct(e2);

    // ray and triangle are parallel if det is close to 0
    det = e1.dotProduct(pvec);
    if (fabs(det) < epsilon) return -1;

    invDet = 1 / det;
    Vec3 tvec = point - v0;
    u = tvec.dotProduct(pvec) * invDet;
    if (u < zero || u > one) return -1;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(e1);
    v = direction.dotProduct(qvec) * invDet;
    if (v < zero || u + v > one) return -1; // Check second & third barycentric coordinate

    // return the distance from point to triangle
    return fabs(e2.dotProduct(qvec) * invDet);
}

// Determine whether the ray intersects with any of the triangles on the surface mesh
bool RaySurfaceIntersect::ray_intersects_surface(const Vec3 &origin, const Vec3 &direction) {
    // Loop through all the faces
    for (int face = 0; face < mesh->faces.size(); ++face)
        if (ray_intersects_triangle(origin, direction, face)) return true;

    return false;
}

// Determine whether the point is whithin cut-off distance from the surface mesh
bool RaySurfaceIntersect::near_surface(const Vec3 &point, const double r_cut) {
    // Loop through all the faces
    for (int face = 0; face < mesh->faces.size(); ++face) {
        double dist = distance_from_triangle(point, face);
        if (dist >= 0 && dist <= r_cut) return true;
        if (dist < 0 && dist >= -0.3*r_cut) return true;
    }

    return false;
}

// Determine whether the point is whithin cut-off distance from the surface mesh
bool RaySurfaceIntersect::near_surface(const Vec3 &point, const Vec3 &direction, const double r_cut) {
    // Loop through all the faces
    for (int face = 0; face < mesh->faces.size(); ++face) {
        double dist = distance_from_triangle(point, direction, face);
        if (dist >= 0 && dist <= r_cut) return true;
    }

    return false;
}

/* ==================================================================
 *  ======================== TetgenMesh ============================
 * ================================================================== */

// Initialize Tetgen data
TetgenMesh::TetgenMesh() {
    tetIOin.initialize();
    tetIOout.initialize();
}

// Group hexahedra around central tetrahedral node
void TetgenMesh::group_hexahedra() {
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // find which hexahedra correspond to which tetrahedral node
    // hexahedra with the same tetrahedral node form the pseudo Voronoi cell of that node
    for (int i = 0; i < hexahedra.size(); ++i) {
        for (int node : hexahedra[i])
            if (node >= node_min && node <= node_max) {
                hexahedra.set_marker(i, node);
                break;
            }
    }
}

// Generate list of nodes that surround the tetrahedral nodes
// The resulting cells resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells
void TetgenMesh::get_pseudo_vorocells(vector<vector<unsigned int>>& cells) const {
    cells = vector<vector<unsigned int>>(nodes.stat.n_tetnode);
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // find the pseudo Voronoi cell nodes for the tetrahedral nodes
    for (int i = 0; i < hexahedra.size(); ++i) {
        const int tetnode = hexahedra.get_marker(i);
        expect(tetnode >= node_min && tetnode <= node_max, "Hexahedron " + to_string(i) +
                " is not marked by the tetrahedral node: " + to_string(tetnode));

        for (unsigned int node : hexahedra[i])
            if ( node != tetnode && nodes.get_marker(node) >= TYPES.EDGECENTROID )
                cells[tetnode].push_back(node);
    }
}

// Function to generate simple mesh that consists of one tetrahedron
bool TetgenMesh::generate_simple() {
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
bool TetgenMesh::recalc() {
    nodes.recalc();
    edges.recalc();
    faces.recalc();
    elems.recalc();
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
bool TetgenMesh::recalc(const string& cmd) {
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
bool TetgenMesh::recalc(const string& cmd1, const string& cmd2) {
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

// Write mesh into files using Tetgen built-in functions
bool TetgenMesh::write(const string& file_name) {
    string file_type = get_file_type(file_name);
    require(file_type == file_name, "Only file names without extension are supported: " + file_type);

    // k - write vtk, Q - quiet, I - suppresses iteration numbers,
    // F - suppress output of .face and .edge, E - suppress output of .ele
    const string cmd = "kIFEQ";
    tetgenbehavior tetgenbeh;

    try {
        tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
        for (unsigned i = 0; i < file_name.size(); ++i)
            tetgenbeh.outfilename[i] = file_name[i];

        tetrahedralize(&tetgenbeh, &tetIOout, NULL);
    } catch (int e) { return 1; }

    return 0;
}

// Function to generate mesh from surface, bulk and vacuum atoms
bool TetgenMesh::generate(const Medium& bulk, const Medium& surf, const Medium& vacuum, const string& cmd) {
    const int n_bulk = bulk.size();
    const int n_surf = surf.size();
    const int n_vacuum = vacuum.size();

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

// Separate tetrahedra into hexahedra
bool TetgenMesh::generate_hexahedra() {
    tethex::Mesh hexmesh;
    hexmesh.read_femocs(this);
    hexmesh.convert();
    hexmesh.export_femocs(this);

    return 0;
}

// Generate manually edges and surface faces
bool TetgenMesh::generate_appendices() {
    // Generate edges from elements
//    generate_edges();
    // Generate surface faces from elements
    generate_surf_faces();
    return 0;
}

// Generate edges from the elements
void TetgenMesh::generate_edges() {
    const int n_elems = elems.size();
    const int n_edges = n_edges_per_elem * n_elems;

    edges.init(n_edges);

    for (SimpleElement selem : elems) {
        for (int e = 0; e < n_edges_per_elem; ++e)
            edges.append(selem.edge(e));
    }
}

// Generate manually surface faces from elements and surface nodes
void TetgenMesh::generate_surf_faces() {
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

//    // Reserve memory for surface faces
//    faces.init( vector_sum(elem_on_surface) + 2 );
//
//    // Make two big faces that pass the xy-min-max corners of surface
//    faces.append(SimpleFace(0, 1, 2));
//    faces.append(SimpleFace(0, 2, 3));

    // Reserve memory for surface faces
    faces.init( vector_sum(elem_on_surface) );

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
bool TetgenMesh::separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd) {
    vector<bool> tet_mask = vector_equal(elems.get_markers(), TYPES.VACUUM);
    vector<bool> hex_mask = vector_equal(hexahedra.get_markers(), TYPES.VACUUM);

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

// Mark mesh nodes and elements
bool TetgenMesh::mark_mesh() {
//    if (mark_nodes())
    if (mark_nodes_vol2())
        return 1;

    mark_elems();
    return 0;
}

// Calculate the neighbourlist for the nodes.
// Two nodes are considered neighbours if they share a tetrahedron.
void TetgenMesh::calc_nborlist(vector<vector<int>>& nborlist) {
    nborlist = vector<vector<int>>(nodes.size());
    for (SimpleElement elem : elems)
        for (int n1 : elem)
            for (int n2 : elem) {
                if (n1 == n2) continue;
                nborlist[n1].push_back(n2);
            }
}

// Mark the nodes by using DBSCAN algorithm (the same as in cluster analysis)
bool TetgenMesh::mark_nodes() {
    int node;
    vector<int> neighbours;
    
    // Calculate neighbourlist for nodes
    vector<vector<int>> nborlist;
    calc_nborlist(nborlist);
            
    // Mark all the nodes with initial values
    nodes.init_markers(nodes.size(), TYPES.NONE);
        
    // Mark the surface, bulk and vacuum nodes by their known position in array
    for (node = nodes.indxs.surf_start; node <= nodes.indxs.surf_end; ++node)
        nodes.set_marker(node, TYPES.SURFACE);
    for (node = nodes.indxs.bulk_start; node <= nodes.indxs.bulk_end; ++node)
        nodes.set_marker(node, TYPES.BULK);
    for (node = nodes.indxs.vacuum_start; node <= nodes.indxs.vacuum_end; ++node)
        nodes.set_marker(node, TYPES.VACUUM);

    // Mark the bulk nodes
    neighbours = nborlist[nodes.indxs.bulk_start];
    for (int i = 0; i < neighbours.size(); ++i) {
        node = neighbours[i];
        if (nodes.get_marker(node) == TYPES.NONE) {
            nodes.set_marker(node, TYPES.BULK);
            neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
        }
    }
    
    // Mark the vacuum nodes
    neighbours = nborlist[nodes.indxs.vacuum_start];
    for (int i = 0; i < neighbours.size(); ++i) {
        node = neighbours[i];
        if (nodes.get_marker(node) == TYPES.NONE) {
            nodes.set_marker(node, TYPES.VACUUM);
            neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
        }
    }
    
    // Nodes inside the thin nanotip may not have nearest neighbour connection
    // with the rest of the bulk. Therefore mark them separately
    for (node = 0; node < nodes.size(); ++node)
        if (nodes.get_marker(node) == TYPES.NONE)
            nodes.set_marker(node, TYPES.BULK);
            
    // Check the result with the number of vacuum atoms: if the surface is too coarse,
    // all the atoms (except the ones added manually to the vacuum)
    // will be marked as either surface or bulk atom
    nodes.calc_statistics();
    expect(nodes.stat.n_vacuum > 4, "Surface is too coarse or rough! Check the output/surface_coarse.xyz,\n"
            "make sure the radius is big enough and consider altering the coarsening factor!");
    return nodes.stat.n_vacuum <= 4;
}

// Calculate factors that show many connections the non-surface nodes have with its non-surface neighbors.
// Nodes with small amount of neighbours are either on the boundary of simubox or on the edge of a hole,
// while nodes with large amount of neighbours are inside the bulk or vacuum domain.
bool TetgenMesh::calc_ranks(vector<int>& ranks, const vector<vector<int>>& nborlist) {
    const int n_nbor_layers = 4;  // number of nearest tetrahedra whose nodes will act as a seed
    const int n_nodes = nodes.size();
    const double max_rank = 100.0;
    double t0;

    // initialise all the ranks to 0
    ranks = vector<int>(n_nodes);

    // distinguish the ranks of surface nodes
    for (int i = nodes.indxs.surf_start; i <= nodes.indxs.surf_end; ++i)
        ranks[i] = -1;

    // calculate the ranks from vacuum side
    vector<int> neighbours = nborlist[nodes.indxs.vacuum_start];
    for (int i = 0; i < neighbours.size(); ++i) {
        int node = neighbours[i];
        if (ranks[node] != -1 && ranks[node]++ == 0)
            neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
    }

    // normalise all the ranks with respect to the maximum rank
    double norm_factor = max_rank / *max_element(ranks.begin(), ranks.end());
    for (int& r : ranks)
        if (r > 0)
            r *= norm_factor;
    
    
    // force the ranks around the vacuum seed region to the maximum value

    vector<vector<int>> nbors(n_nbor_layers);
    ranks[nodes.indxs.vacuum_start] = max_rank;
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of node neighbour list
        if (layer == 0)
            nbors[0] = nborlist[nodes.indxs.vacuum_start];
        else {
            for (int nbor : nbors[layer-1])
                if (ranks[nbor] > 0)
                    nbors[layer].insert(nbors[layer].end(), nborlist[nbor].begin(), nborlist[nbor].end());
        }
        for (int nbor : nbors[layer])
            if (ranks[nbor] > 0)
                ranks[nbor] = max_rank;
    }

    // check whether only the vacuum nodes were ranked or did the ranks penetrate trough the surface
    // the latter suggests, that the tetrahedral surface has holes inside that may or may not cause problems
    for (int r : ranks)
        if (r == 0)
            return true;

    return false;
}

// Mark the nodes by using the DBSCAN clustering algorithm
// To increase the tolerance against the holes in the surface, calculate the ranking of the nodes,
// i.e the measure how close the node is to the hole or to the simubox boundary.
// Nodes with small ranking act as a boundary for the clustering algorithm.
// The ranking helps to get rid of the effect of small holes.
bool TetgenMesh::mark_nodes_vol2() {
    const int n_nodes = nodes.size();
    const int min_rank = 33;
    int node;

    // Calculate neighbour list for nodes
    vector<vector<int>> nborlist;
    calc_nborlist(nborlist);

    // Calculate the ranks for the nodes to increase the robustness of the bulk-vacuum separator
    vector<int> ranks;
    if ( !calc_ranks(ranks, nborlist) )
        write_silent_msg("Surface has holes, therefore the mesh may or might not be valid!\n"
                "Check the output/hexmesh_bulk.vtk and in case of problems, make sure \n"
                "the radius is big enough and consider altering the coarsening factors!");

    // Mark all the nodes with initial values
    nodes.init_markers(n_nodes, TYPES.NONE);

    // Mark the surface, bulk and vacuum nodes by their known position in array
    for (node = nodes.indxs.surf_start; node <= nodes.indxs.surf_end; ++node)
        nodes.set_marker(node, TYPES.SURFACE);

    for (node = nodes.indxs.bulk_start; node <= nodes.indxs.bulk_end; ++node)
        nodes.set_marker(node, TYPES.BULK);

    for (node = nodes.indxs.vacuum_start; node <= nodes.indxs.vacuum_end; ++node)
        nodes.set_marker(node, TYPES.VACUUM);

    // Mark the vacuum nodes
    vector<int> neighbours = nborlist[nodes.indxs.vacuum_start];
    for (int i = 0; i < neighbours.size(); ++i) {
        node = neighbours[i];
        if (nodes.get_marker(node) == TYPES.NONE) {
            nodes.set_marker(node, TYPES.VACUUM);
            if (ranks[node] >= min_rank)
                neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
        }
    }

    // Mark the bulk nodes
    // no need to check the node ranks as vacuum nodes are all already marked
    neighbours = nborlist[nodes.indxs.bulk_start];
    for (int i = 0; i < neighbours.size(); ++i) {
        node = neighbours[i];
        if (nodes.get_marker(node) == TYPES.NONE) {
            nodes.set_marker(node, TYPES.BULK);
            neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
        }
    }

    // Nodes inside the thin nanotip may not have nearest neighbour connection
    // with the rest of the bulk. Therefore mark them separately
    for (int &marker : *nodes.get_markers())
        if (marker == TYPES.NONE)
            marker = TYPES.BULK;

    return 0;
}

// Locate the tetrahedron by the location of its nodes
int TetgenMesh::locate_element(SimpleElement& elem) {
    // Element in vacuum is supposed not to have nodes in bulk,
    //  in bulk is supposed not to have nodes in vacuum,
    //  in surface is supposed to consist only of nodes in surface
    for (int node : elem) {
        const int nodemarker = nodes.get_marker(node);
        if (nodemarker == TYPES.VACUUM)
            return TYPES.VACUUM;
        else if (nodemarker == TYPES.BULK)
            return TYPES.BULK;
    }
    
    return TYPES.SURFACE;
}

// Mark the tetrahedra by the location of its nodes
void TetgenMesh::mark_elems() {
    // Reserve memory for markers
    elems.init_markers(elems.size());

    // Locate all the elements
    for (SimpleElement elem : elems)
        elems.append_marker(locate_element(elem));
}

// Mark the edges on the simulation cell perimeter by the node markers
void TetgenMesh::mark_edges() {
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
void TetgenMesh::mark_faces() {
    const double eps = 0.1 * elems.stat.edgemin;
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
