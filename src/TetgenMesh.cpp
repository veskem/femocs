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

// Initialize Tetgen data
TetgenMesh::TetgenMesh() {
    tetIOin.initialize();
    tetIOout.initialize();
}

void TetgenMesh::test_mapping() const {
    cout << "\ntets of tris:" << endl;
    for (int i = 0; i < 10; ++i) {
        cout << i << ":\t";
        for (int tet : faces.to_tets(i))
            cout << tet << ", ";
        cout << endl;
    }

    cout << "tris of tets:" << endl;
    for (int i = 0; i < 10; ++i) {
        cout << i << ":\t";
        for (int tri : elems.to_tris(i))
            cout << tri << ", ";
        cout << endl;
    }

    cout << "hexs of quads:" << endl;
    for (int i = 0; i < 10; ++i) {
        cout << i << ":\t";
        for (int hex : quads.to_hexs(i))
            cout << hex << ", ";
        cout << endl;
    }

    cout << "quads of hexs:" << endl;
    for (int i = 0; i < 10; ++i) {
        cout << i << ":\t";
        for (int quad : hexahedra.to_quads(i))
            cout << quad << ", ";
        cout << endl;
    }

    int n_cells, n_mapped_cells;

    n_cells = faces.size();
    n_mapped_cells = 0;
    for (int i = 0; i < n_cells; ++i)
        n_mapped_cells += faces.to_tets(i).size();
    printf("faces: n_cells=%i, n_mapped_cells=%i\n", n_cells, n_mapped_cells);

    n_cells = elems.size();
    n_mapped_cells = 0;
    for (int i = 0; i < n_cells; ++i)
        n_mapped_cells += elems.to_tris(i).size();
    printf("elems: n_cells=%i, n_mapped_cells=%i\n", n_cells, n_mapped_cells);

    n_cells = quads.size();
    n_mapped_cells = 0;
    for (int i = 0; i < n_cells; ++i)
        n_mapped_cells += quads.to_hexs(i).size();
    printf("quads: n_cells=%i, n_mapped_cells=%i\n", n_cells, n_mapped_cells);

    n_cells = hexahedra.size();
    n_mapped_cells = 0;
    for (int i = 0; i < n_cells; ++i)
        n_mapped_cells += hexahedra.to_quads(i).size();
    printf("hexahedra: n_cells=%i, n_mapped_cells=%i\n", n_cells, n_mapped_cells);
}

int TetgenMesh::tri2tet(const int tri, const int region) const {

    if (region == TYPES.VACUUM) {
        for (int tet : faces.to_tets(tri))
            if (tet >= 0 && elems.get_marker(tet) == TYPES.VACUUM)
                return tet;
    } else if (region == TYPES.BULK) {
        for (int tet : faces.to_tets(tri))
            if (tet >= 0 && elems.get_marker(tet) != TYPES.VACUUM)
                return tet;
    } else
        require(false, "Unimplemented region: " + to_string(region));

    return -1;
}

// map quadrangle to hexahedron that is located in a given region
int TetgenMesh::quad2hex(const int quad, const int region) const {

    if (region == TYPES.VACUUM) {
        for (int hex : quads.to_hexs(quad))
            if (hex >= 0 && elems.get_marker(hexahedra.to_tet(hex)) == TYPES.VACUUM)
                return hex;
    } else if (region == TYPES.BULK) {
        for (int hex : quads.to_hexs(quad))
            if (hex >= 0 && elems.get_marker(hexahedra.to_tet(hex)) != TYPES.VACUUM)
                return hex;
    } else
        require(false, "Unimplemented region: " + to_string(region));

    return -1;
}

// Delete the data of previously stored mesh and initialise a new one
void TetgenMesh::clear() {
    tetIOin.deinitialize();
    tetIOout.deinitialize();
    tetIOin.initialize();
    tetIOout.initialize();
}

// Smoothen the triangles using different versions of Taubin smoothing algorithm
// Code is inspired from the work of Shawn Halayka
// TODO if found, add the link to Shawn's version
void TetgenMesh::smoothen(const int n_steps, const double lambda, const double mu, const string& algorithm) {
    if (algorithm == "none") return;

    // determine the neighbouring between nodes connected to the faces
    vector<vector<unsigned>> nborlist;
    faces.calc_nborlist(nborlist);

    // remove the nodes that are on the boundary of simubox
    const double eps = 0.01 * elems.stat.edgemin;
    nodes.calc_statistics();
    for (int i = 0; i < nodes.size(); ++i) {
        Point3 p = nodes[i];
        if (on_boundary(p.x, nodes.stat.xmin, nodes.stat.xmax, eps) ||
                on_boundary(p.y, nodes.stat.ymin, nodes.stat.ymax, eps) ||
                on_boundary(p.z, nodes.stat.zmin, nodes.stat.zmax, eps))
            nborlist[i] = vector<unsigned>();
    }

    // run the Taubin smoothing
    if (algorithm == "laplace") {
        for (size_t s = 0; s < n_steps; ++s) {
            laplace_smooth(lambda, nborlist);
            laplace_smooth(mu, nborlist);
        }
    }

    else if (algorithm == "fujiwara") {
        for (size_t s = 0; s < n_steps; ++s) {
            fujiwara_smooth(lambda, nborlist);
            fujiwara_smooth(mu, nborlist);
        }
    }

    faces.calc_norms_and_areas();
}

// Smoothen the surface mesh using Taubin lambda|mu algorithm with inverse neighbour count weighting
void TetgenMesh::laplace_smooth(const double scale, const vector<vector<unsigned>>& nborlist) {
    size_t n_nodes = nodes.size();
    vector<Point3> displacements(n_nodes);

    // Get per-vertex displacement
    for(size_t i = 0; i < n_nodes; ++i) {
        // Skip vertices that are not on the surface
        if (nborlist[i].size() == 0)
            continue;

        const double weight = 1.0 / nborlist[i].size();

        // Sum the displacements
        for(size_t nbor : nborlist[i])
            displacements[i] += (nodes[nbor] - nodes[i]) * weight;
    }

    // Apply per-vertex displacement
    for (size_t i = 0; i < n_nodes; ++i)
        nodes.set_node(i, nodes[i] + displacements[i]*scale);
}

// Smoothen the surface mesh using Taubin lambda|mu algorithm with Fujiwara weighting
void TetgenMesh::fujiwara_smooth(const double scale, const vector<vector<unsigned>>& nborlist) {
    vector<Point3> displacements(nodes.size());

    // Get per-vertex displacement.
    for (size_t i = 0; i < nodes.size(); ++i) {
        // Skip vertices that are not on the surface.
        if (nborlist[i].size() == 0)
            continue;

        vector<double> weights(nborlist[i].size());
        Point3 node = nodes[i];

        // Calculate Fujiwara weights based on edge lengths.
        for (size_t j = 0; j < nborlist[i].size(); j++) {
            size_t nbor = nborlist[i][j];

            double edge_length = node.distance(nodes[nbor]);
            expect(edge_length > 0, "Zero edge not allowed!");

            weights[j] = 1.0 / edge_length;
        }

        // Normalize the weights so that they sum up to 1.
        double s = vector_sum(weights);
        if (s == 0) s = numeric_limits<float>::epsilon();

        s = 1.0 / s;
        for (size_t j = 0; j < weights.size(); j++)
            weights[j] *= s;

        // Sum the displacements.
        for(size_t j = 0; j < nborlist[i].size(); j++) {
            size_t nbor = nborlist[i][j];
            displacements[i] += (nodes[nbor] - node) * weights[j];
        }
    }

    // Apply per-vertex displacement.
    for (size_t i = 0; i < nodes.size(); i++)
        nodes.set_node(i, nodes[i] + displacements[i]*scale);
}

// Smoothen the surface mesh using Taubin lambda|mu algorithm with curvature normal weighting
void TetgenMesh::curvature_norm_smooth(const double scale, const vector<vector<unsigned>>& nborlist) {
    size_t n_nodes = nodes.size(); 
    vector<Point3> displacements(n_nodes);

    // Generate vertex to triangle mapping
    vector<vector<int>> vertex_to_triangle_indices = vector<vector<int>>(n_nodes);
    int f = 0;
    for (SimpleFace face : faces)
        for (int node : face)
            vertex_to_triangle_indices[node].push_back(f++);

    // Get per-vertex displacement
    for(size_t i = 0; i < n_nodes; ++i) {
        if (nborlist[i].size() == 0)
            continue;

        vector<double> weights(nborlist[i].size());
        size_t angle_error = 0;

        // For each vertex pair (ie. each edge),
        // calculate weight based on the two opposing angles (ie. curvature normal scheme).
        for (size_t j = 0; j < nborlist[i].size(); ++j) {
            size_t nbor = nborlist[i][j];
            size_t angle_count = 0;

            // Find out which two triangles are shared by the edge.
            for (size_t tri0_index : vertex_to_triangle_indices[i]) {
                for (size_t tri1_index : vertex_to_triangle_indices[nbor]) {

                    // tri0_index == tri1_index will occur twice per edge.
                    if (tri0_index != tri1_index) continue;

                    // Find the third vertex in this triangle (the vertex that doesn't belong to the edge).
                    for (size_t opp_vert_index : faces[tri0_index]) {
                        // This will occur once per triangle.
                        if (opp_vert_index != i && opp_vert_index != nbor) {
                            // Get the angle opposite of the edge.
                            Vec3 c = nodes.get_vec(opp_vert_index);
                            Vec3 a = nodes.get_vec(i)    - c;
                            Vec3 b = nodes.get_vec(nbor) - c;
                            a.normalize();
                            b.normalize();

                            double dot = a.dotProduct(b);
                            dot = min(1.0, max(-1.0, dot));
                            double angle = acos(dot);

                            // Curvature normal weighting.
                            double slope = tan(angle);

                            if(slope == 0)
                                slope = numeric_limits<double>::epsilon();

                            // Note: Some weights will be negative, due to obtuse triangles.
                            // You may wish to do weights[j] += fabsf(1.0f / slope); here.
                            weights[j] += 1.0 / slope;

                            angle_count++;

                            break;
                        }
                    }

                    // Since we found a triangle match, we can skip to the first vertex's next triangle.
                    break;
                }
            } // End of: Find out which two triangles are shared by the vertex pair.

            if (angle_count != 2)
                angle_error++;

        } // End of: For each vertex pair (ie. each edge).

//        if (angle_error != 0) {
//            ostringstream os;
//            os<< "Warning: Vertex " << i << " belongs to " << angle_error
//                    << " edges that do not belong to two triangles ("
//                    << nborlist[i].size() - angle_error << " edges were OK).\n"
//                    << "Your mesh probably has cracks or holes in it." << endl;
//            write_silent_msg(os.str());
//        }

        // Normalize the weights so that they sum up to 1.

        // Note: Some weights will be negative, due to obtuse triangles.
        // You may wish to add together absolute values of weights.
        double s = 0;
        for (double w : weights)
            s += fabs(w);

        if (s == 0) s = numeric_limits<float>::epsilon();

        s = 1.0 / s;
        for (size_t j = 0; j < weights.size(); j++)
            weights[j] *= s;

        // Sum the displacements.
        for (size_t j = 0; j < nborlist[i].size(); j++) {
            size_t nbor = nborlist[i][j];
            displacements[i] += (nodes[nbor] - nodes[i]) * weights[j];
        }
    }

    // TODO: Find out why there are cases where displacement is much, much, much larger than all edge lengths put together.

    // Apply per-vertex displacement.
    for (size_t i = 0; i < n_nodes; i++)
        nodes.set_node(i, nodes[i] + displacements[i] * scale);
}

// Function to generate simple mesh that consists of one tetrahedron
int TetgenMesh::generate_simple() {
    const int n_nodes = elems.DIM;
    const int n_elems = 1;

    nodes.init(n_nodes);
    nodes.append(Point3(1.0, 0.0, 0.7));
    nodes.append(Point3(-1.0, 0.0, 0.7));
    nodes.append(Point3(0.0, 1.0, -0.7));
    nodes.append(Point3(0.0, -1.0, -0.7));

    faces.init(n_tris_per_tet);
    faces.append(SimpleFace(0, 1, 3));
    faces.append(SimpleFace(1, 2, 3));
    faces.append(SimpleFace(2, 0, 3));
    faces.append(SimpleFace(0, 1, 2));

    edges.init(n_edges_per_tet);
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

// Copy mesh from input to output or vice versa without modification
int TetgenMesh::transfer(const bool write2read) {
    nodes.transfer(write2read);
    edges.transfer(write2read);
    faces.transfer(write2read);
    elems.transfer(write2read);
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
int TetgenMesh::recalc(const string& cmd) {
    try {
        tetrahedralize(const_cast<char*>(cmd.c_str()), &tetIOin, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        faces.set_counter(tetIOout.numberoftrifaces);
        elems.set_counter(tetIOout.numberoftetrahedra);
        faces.calc_statistics();
        elems.calc_statistics();
    }
    catch (int e) { return e; }
    return 0;
}

// Function to perform tetgen calculation on input and write_tetgen data
int TetgenMesh::recalc(const string& cmd1, const string& cmd2) {
    try {
        tetgenio tetIOtemp;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp);
        tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        faces.set_counter(tetIOout.numberoftrifaces);
        elems.set_counter(tetIOout.numberoftetrahedra);
        faces.calc_statistics();
        elems.calc_statistics();
    }
    catch (int e) { return e; }
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
int TetgenMesh::generate(const Medium& bulk, const Medium& surf, const Medium& vacuum, const string& cmd) {
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

// Group hexahedra & quadrangles around central tetrahedral & triangular nodes
void TetgenMesh::group_hexahedra() {
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;
    const int n_hexs = hexahedra.size();
    const int n_quads = quads.size();

    // find which hexahedra correspond to which tetrahedral node
    // hexahedra with the same tetrahedral node form the pseudo 3D Voronoi cell of that node
    for (int i = 0; i < hexahedra.size(); ++i)
        // group hexs in vacuum
        if (hexahedra.get_marker(i) == TYPES.VACUUM) {
            for (int node : hexahedra[i])
                if (node >= node_min && node <= node_max) {
                    hexahedra.set_marker(i, 1 + node);
                    break;
                }
        // group hexs in bulk
        } else {
            for (int node : hexahedra[i])
                if (node >= node_min && node <= node_max) {
                    hexahedra.set_marker(i, -1 - node);
                    break;
                }
        }

    // find which quadrangle correspond to which triangular node
    // quadrangles with the same triangular node form the pseudo 2D Voronoi cell of that node
    for (int i = 0; i < quads.size(); ++i) {
        for (int node : quads[i])
            if (node >= node_min && node <= node_max) {
                quads.set_marker(i, node);
                break;
            }
    }
}

// Separate tetrahedra & triangles into hexahedra & quadrangles
bool TetgenMesh::generate_hexahedra() {
    tethex::Mesh hexmesh;
    hexmesh.read_femocs(this);
    hexmesh.convert();
    hexmesh.export_femocs(this);

    group_hexahedra();
    calc_quad2hex_mapping();
    nodes.calc_statistics();

    return 0;
}

// Using the separated tetrahedra generate the triangular surface on the vacuum-material boundary
int TetgenMesh::generate_surface(const Medium::Sizes& sizes, const string& cmd1, const string& cmd2) {
    TetgenMesh vacuum;
    vector<bool> tet_mask = vector_equal(elems.get_markers(), TYPES.VACUUM);

    // Transfer vacuum nodes and tetrahedra
    vacuum.nodes.copy(this->nodes);
    vacuum.elems.copy(this->elems, tet_mask);

    // calculate surface triangles
    int error_code = vacuum.recalc(cmd1);
    if (error_code) return error_code;

    // copy the triangles that are not on the simubox perimeter to input tetIO
    vacuum.faces.calc_statistics();
    const int n_surf_faces = faces.copy_surface(vacuum.faces, sizes);

    // transfer elements and nodes to input
    nodes.transfer(false);
    elems.transfer(false);

    // calculate the tetrahedron-triangle connectivity
    // the simubox boundary faces must also be calculated, no way to opt-out
    error_code = recalc(cmd2);
    if (error_code) return error_code;

    // the faces on the simubox sides are appended after the surface faces.
    // such property allows to remove the faces on the sides without affecting the tri2tet mapping.
    // such cleaning is useful to make other workflow faster.
    faces.init(n_surf_faces);
    for (int i = 0; i < n_surf_faces; ++i)
        faces.append(faces[i]);
    faces.transfer();

    calc_tet2tri_mapping();
    return 0;
}

void TetgenMesh::calc_tet2tri_mapping() {
    const int n_tris = faces.size();
    vector<vector<int>> tet2tri_map(elems.size());

    for (int i = 0; i < n_tris; ++i) {
        for (int tet : faces.to_tets(i))
            if (tet >= 0)
                tet2tri_map[tet].push_back(i);
    }
    elems.store_map(tet2tri_map);
}

void TetgenMesh::calc_quad2hex_mapping() {
    const int n_quads = quads.size();
    const int n_hexs = hexahedra.size();

    vector<array<int,2>> quad2hex_map = vector<array<int,2>>(n_quads, {-1,-1});
    vector<vector<int>> hex2quad_map = vector<vector<int>>(n_hexs);

    for (int quad = 0; quad < n_quads; ++quad) {
        SimpleQuad squad = quads[quad];

        // loop through the tetrahedra that are connected to the quadrangle
        int region = 0;
        for (int tet : faces.to_tets(quads.to_tri(quad))) {
            if (tet < 0) continue;

            // loop through all the hexahedra connected to the tetrahedron
            for (int hex : elems.to_hexs(tet)) {

                // count for the # common nodes between quadrangle and hexahedron
                int n_common_nodes = 0;
                for (unsigned int node : hexahedra[hex])
                    n_common_nodes += squad == node;

                // quad belongs to hex, if they share 4 nodes
                if (n_common_nodes == n_nodes_per_quad) {
                    quad2hex_map[quad][region++] = hex;
                    hex2quad_map[hex].push_back(quad);
                }
            }
        }
    }

    // store the mapping on the cells side
    quads.store_map(quad2hex_map);
    hexahedra.store_map(hex2quad_map);
}

// Generate manually surface faces from elements and surface nodes
void TetgenMesh::generate_manual_surface() {
    const int n_elems = elems.size();
    const int max_surf_indx = nodes.indxs.surf_end;

    // booleans showing whether element i has exactly one face on the surface or not
    vector<bool> elem_on_surface; elem_on_surface.reserve(n_elems);
    // booleans showing whether node i is on the surface or not
    vector<bool> surf_locs(4);

    // Mark the elements that have exactly one face on the surface
    for (SimpleElement elem : elems) {
        for (int i = 0; i < 4; ++i)
            surf_locs[i] = elem[i] <= max_surf_indx;
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
            for (int i = 0; i < 4; ++i)
                surf_locs[i] = elem[i] <= max_surf_indx;

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

// Generate edges from the elements
void TetgenMesh::generate_edges() {
    const int n_elems = elems.size();
    edges.init(n_edges_per_tet * n_elems);

    for (SimpleElement selem : elems) {
        for (int e = 0; e < n_edges_per_tet; ++e)
            edges.append(selem.edge(e));
    }
}

// Separate vacuum and bulk mesh from the union mesh by the element markers
int TetgenMesh::separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd) {
    vector<bool> tet_mask = vector_equal(elems.get_markers(), TYPES.VACUUM);
    vector<bool> hex_mask = vector_equal(hexahedra.get_markers(), TYPES.VACUUM);

    // Transfer vacuum nodes, tetrahedra, hexahedra and their markers
    vacuum.nodes.copy(this->nodes);
    vacuum.nodes.copy_markers(this->nodes);
    vacuum.faces.copy(this->faces);
    vacuum.faces.copy_markers(this->faces);
    vacuum.quads.copy(this->quads);
    vacuum.quads.copy_markers(this->quads);
    vacuum.elems.copy(this->elems, tet_mask);
    vacuum.elems.copy_markers(this->elems, tet_mask);
    vacuum.hexahedra.copy(this->hexahedra, hex_mask);
    vacuum.hexahedra.copy_markers(this->hexahedra, hex_mask);

    tet_mask.flip();
    hex_mask.flip();

    // Transfer bulk nodes, tetrahedra, hexahedra and their markers
    bulk.nodes.copy(this->nodes);
    bulk.nodes.copy_markers(this->nodes);
    bulk.faces.copy(this->faces);
    bulk.faces.copy_markers(this->faces);
    bulk.quads.copy(this->quads);
    bulk.quads.copy_markers(this->quads);
    bulk.elems.copy(this->elems, tet_mask);
    bulk.elems.copy_markers(this->elems, tet_mask);
    bulk.hexahedra.copy(this->hexahedra, hex_mask);
    bulk.hexahedra.copy_markers(this->hexahedra, hex_mask);

    return vacuum.recalc(cmd) + bulk.recalc(cmd);
}

// Write bulk or vacuum mesh
void TetgenMesh::write_separate(const string& file_name, const int type) {
    vector<bool> hex_mask;
    if (type == TYPES.VACUUM)
        hex_mask = vector_greater(hexahedra.get_markers(), 0);
    else
        hex_mask = vector_less(hexahedra.get_markers(), 0);

    TetgenMesh tempmesh;
    tempmesh.nodes.copy(this->nodes);
    tempmesh.nodes.copy_markers(this->nodes);
    tempmesh.nodes.transfer();
    tempmesh.hexahedra.copy(this->hexahedra, hex_mask);
    tempmesh.hexahedra.copy_markers(this->hexahedra, hex_mask);

    tempmesh.hexahedra.write(file_name);
}

// Mark mesh nodes and elements
bool TetgenMesh::mark_mesh() {
//    if (mark_nodes())
    if (rank_and_mark_nodes())
        return 1;

    mark_elems();
    return 0;
}

// Generate list of nodes that surround the tetrahedral nodes
// The resulting cells resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells
void TetgenMesh::calc_pseudo_3D_vorocells(vector<vector<unsigned>>& cells, const bool vacuum) const {
    cells = vector<vector<unsigned>>(nodes.stat.n_tetnode);
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // the sign of hexahedron marker shows its location (>0 == vacuum, <0 == bulk)
    // and the (magnitude - 1) the index of tetrahedral node it is connected to
    int multiplier = 1;
    if (!vacuum) multiplier = -1;

    // find the pseudo Voronoi cell nodes for the tetrahedral nodes
    for (int hex = 0; hex < hexahedra.size(); ++hex) {
        int tetnode = multiplier * hexahedra.get_marker(hex);
        if (tetnode <= 0) continue;
        tetnode -= 1;

        for (int node : hexahedra[hex])
            if ( node != tetnode && nodes.get_marker(node) >= TYPES.EDGECENTROID )
                cells[tetnode].push_back(node);
    }
}

// Generate list of quadrangle nodes that surround the triangle nodes
// The resulting cells resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells
void TetgenMesh::calc_pseudo_2D_vorocells(vector<vector<unsigned>>& cells) const {
    cells = vector<vector<unsigned>>(nodes.stat.n_tetnode);
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    for (int quad = 0; quad < quads.size(); ++quad) {
        const int trinode = quads.get_marker(quad);
        expect(trinode >= node_min && trinode <= node_max, "Quadrangle " + to_string(quad) +
                " is not marked by the triangle node: " + to_string(trinode));

        for (int node : quads[quad])
            if (node != trinode)
                cells[trinode].push_back(node);
    }
}

// Mark the nodes by using DBSCAN algorithm (the same as in cluster analysis)
bool TetgenMesh::mark_nodes() {
    int node;
    vector<unsigned> neighbours;
    
    // Calculate neighbourlist for nodes
    vector<vector<unsigned>> nborlist;
    elems.calc_nborlist(nborlist);
            
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
    for (size_t i = 0; i < neighbours.size(); ++i) {
        node = neighbours[i];
        if (nodes.get_marker(node) == TYPES.NONE) {
            nodes.set_marker(node, TYPES.BULK);
            neighbours.insert(neighbours.end(), nborlist[node].begin(), nborlist[node].end());
        }
    }
    
    // Mark the vacuum nodes
    neighbours = nborlist[nodes.indxs.vacuum_start];
    for (size_t i = 0; i < neighbours.size(); ++i) {
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
    expect(nodes.stat.n_vacuum > 4, "Surface is too coarse or rough! Check the out/surface_coarse.xyz,\n"
            "make sure the radius is big enough and consider altering the coarsening factor!");
    return nodes.stat.n_vacuum <= 4;
}

// Calculate factors that show many connections the non-surface nodes have with its non-surface neighbors.
// Nodes with small amount of neighbours are either on the boundary of simubox or on the edge of a hole,
// while nodes with large amount of neighbours are inside the bulk or vacuum domain.
bool TetgenMesh::calc_ranks(vector<int>& ranks, const vector<vector<unsigned>>& nborlist) {
    const int n_nbor_layers = 4;  // number of nearest tetrahedra whose nodes will act as a seed
    const int n_nodes = nodes.size();
    const double max_rank = 100.0;
    const double eps = 0.01 * elems.stat.edgemin;

    // initialise all the ranks to 0
    ranks = vector<int>(n_nodes);

    // distinguish the ranks of surface nodes
    for (int i = nodes.indxs.surf_start; i <= nodes.indxs.surf_end; ++i)
        ranks[i] = -1;

    // calculate the ranks from vacuum side
    vector<unsigned> neighbours = nborlist[nodes.indxs.vacuum_start];
    for (size_t i = 0; i < neighbours.size(); ++i) {
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

    vector<vector<unsigned>> nbors(n_nbor_layers);
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

    // force the ranks on the simubox vacuum perimeter to maximum value

    nodes.calc_statistics();
    for (int i = 0; i < n_nodes; ++i) {
        if ( ranks[i] > 0 && (on_boundary(nodes[i].x, nodes.stat.xmin, nodes.stat.xmax, eps) ||
                on_boundary(nodes[i].y, nodes.stat.ymin, nodes.stat.ymax, eps)) )
            ranks[i] = max_rank;
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
bool TetgenMesh::rank_and_mark_nodes() {
    const int n_nodes = nodes.size();
    const int min_rank = 30;
    int node;

    // Calculate neighbour list for nodes
    vector<vector<unsigned>> nborlist;
    elems.calc_nborlist(nborlist);

    // Calculate the ranks for the nodes to increase the robustness of the bulk-vacuum separator
    vector<int> ranks;
    if ( !calc_ranks(ranks, nborlist) )
        write_silent_msg("Surface has holes, therefore the mesh may or might not be valid!\n"
                "Check the out/hexmesh_bulk.vtk and in case of problems, make sure \n"
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
    vector<unsigned> neighbours = nborlist[nodes.indxs.vacuum_start];
    for (size_t i = 0; i < neighbours.size(); ++i) {
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
    for (size_t i = 0; i < neighbours.size(); ++i) {
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
