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

void TetgenMesh::set_write_period(const double dt) {
    FileWriter::set_write_period(dt);
    nodes.set_write_period(dt);
    edges.set_write_period(dt);
    tris.set_write_period(dt);
    tets.set_write_period(dt);
    quads.set_write_period(dt);
    hexs.set_write_period(dt);
}

// Code is inspired from the work of Shawn Halayka
// TODO if found, add the link to Shawn's version
void TetgenMesh::smoothen(const int n_steps, const double lambda, const double mu, const string& algorithm) {
    if (algorithm == "none") return;

    // determine the neighbouring between nodes connected to the faces
    vector<vector<unsigned>> nborlist;
    tris.calc_nborlist(nborlist);

    // remove the nodes that are on the boundary of simubox
    const double eps = 0.01 * tets.stat.edgemin;
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
        for (int s = 0; s < n_steps; ++s) {
            laplace_smooth(lambda, nborlist);
            laplace_smooth(mu, nborlist);
        }
    }

    else if (algorithm == "fujiwara") {
        for (int s = 0; s < n_steps; ++s) {
            fujiwara_smooth(lambda, nborlist);
            fujiwara_smooth(mu, nborlist);
        }
    }
}

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

void TetgenMesh::fujiwara_smooth(const double scale, const vector<vector<unsigned>>& nborlist) {
    vector<Point3> displacements(nodes.size());
    const unsigned int n_nodes = nodes.size();

    // Get per-vertex displacement.
    for (size_t i = 0; i < n_nodes; ++i) {
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
    for (size_t i = 0; i < n_nodes; i++)
        nodes.set_node(i, nodes[i] + displacements[i]*scale);
}

int TetgenMesh::generate_simple() {
    const int n_nodes = tets.DIM;
    const int n_elems = 1;

    nodes.init(n_nodes);
    nodes.append(Point3(1.0, 0.0, 0.7));
    nodes.append(Point3(-1.0, 0.0, 0.7));
    nodes.append(Point3(0.0, 1.0, -0.7));
    nodes.append(Point3(0.0, -1.0, -0.7));

    tris.init(n_tris_per_tet);
    tris.append(SimpleFace(0, 1, 3));
    tris.append(SimpleFace(1, 2, 3));
    tris.append(SimpleFace(2, 0, 3));
    tris.append(SimpleFace(0, 1, 2));

    edges.init(n_edges_per_tet);
    edges.append(SimpleEdge(0, 1));
    edges.append(SimpleEdge(0, 2));
    edges.append(SimpleEdge(0, 3));
    edges.append(SimpleEdge(1, 2));
    edges.append(SimpleEdge(1, 3));
    edges.append(SimpleEdge(2, 3));

    tets.init(n_elems);
    tets.append(SimpleElement(0, 1, 2, 3));

    return recalc("rQ");
}

int TetgenMesh::generate_union(const Medium& bulk, const Medium& surf, const Medium& vacuum, const string& cmd) {
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

int TetgenMesh::generate(const Medium& bulk, const Medium& surf, const Medium& vacuum, const Config& conf) {
    require(bulk.size() > 0,   "Empty mesh generators in bulk detected!");
    require(surf.size() > 0,   "Empty mesh generators on surface detected!");
    require(vacuum.size() > 0, "Empty mesh generators in vacuum detected!");

    // Delete previous mesh
    clear();

    // r - reconstruct, n(n) - output tet neighbour list (and tri-tet connection),
    // Q - quiet, q - mesh quality, a - element volume, E - suppress output of elements
    // F - suppress output of faces and edges, B - suppress output of boundary info
    string command = "rQFBq" + conf.geometry.mesh_quality;
    if (conf.geometry.element_volume != "") command += "a" + conf.geometry.element_volume;

    // Make union mesh where both vacuum and material domain are present
    int err_code = generate_union(bulk, surf, vacuum, command);
    check_return(err_code, "Tetrahedrization failed with error code " + d2s(err_code));

    // Mark tetrahedral mesh
    bool fail = mark_mesh();
    check_return(fail, "Mesh marking failed!");

    // Generate surface faces
    err_code = generate_surface("rQB", "rQnn");
    check_return(err_code, "Triangulation failed with error code " + d2s(err_code));

    // Smoothen surface faces
    if (conf.smoothing.algorithm != "none" && conf.smoothing.n_steps > 0)
        smoothen(conf.smoothing.n_steps, conf.smoothing.lambda_mesh, conf.smoothing.mu_mesh, conf.smoothing.algorithm);

    // has to be separate to ensure that data is for sure calculated
    tris.calc_appendices();

    // Convert tetrahedra & triangles to hexahedra & quadrangles
    fail = generate_hexahedra();
    return fail;
}

int TetgenMesh::read(const string &file_name, const string &cmd) {

    // delete available mesh data
    clear();

    tethex::Mesh hexmesh;
    // read tetrahedral mesh from file
    hexmesh.read(file_name, 0);
    // generate hexahedra and quadrangles
    hexmesh.convert();
    // export mesh to Femocs
    hexmesh.export_all_mesh(this, true);

    // calculate mapping between triangles and tetrahedra
    int err_code = calc_tet2tri_mapping(cmd, tris.get_n_markers());
    check_return(err_code, "Generation of tri2tet mapping failed with error code " + d2s(err_code));

    // calculate mapping between quadrangles and hexahedra
    group_hexahedra();
    calc_quad2hex2quad_mapping();

    nodes.calc_statistics();
    tris.calc_appendices();

    return 0;
}

void TetgenMesh::clear() {
    tetIOin.deinitialize();
    tetIOout.deinitialize();
    tetIOin.initialize();
    tetIOout.initialize();
}

int TetgenMesh::transfer(const bool write2read) {
    nodes.transfer(write2read);
    edges.transfer(write2read);
    tris.transfer(write2read);
    tets.transfer(write2read);
    return 0;
}

int TetgenMesh::recalc(const string& cmd) {
    try {
        tetrahedralize(const_cast<char*>(cmd.c_str()), &tetIOin, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        tris.set_counter(tetIOout.numberoftrifaces);
        tets.set_counter(tetIOout.numberoftetrahedra);
        tris.calc_statistics();
        tets.calc_statistics();
    }
    catch (int e) { return e; }
    return 0;
}

int TetgenMesh::recalc(const string& cmd1, const string& cmd2) {
    try {
        tetgenio tetIOtemp;
        tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp);
        tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp, &tetIOout);
        nodes.set_counter(tetIOout.numberofpoints);
        edges.set_counter(tetIOout.numberofedges);
        tris.set_counter(tetIOout.numberoftrifaces);
        tets.set_counter(tetIOout.numberoftetrahedra);
        tris.calc_statistics();
        tets.calc_statistics();
    }
    catch (int e) { return e; }
    return 0;
}

void TetgenMesh::group_hexahedra() {
    const int node_min = nodes.indxs.tetnode_start;
    const int node_max = nodes.indxs.tetnode_end;

    // find which hexahedra correspond to which tetrahedral node
    // hexahedra with the same tetrahedral node form the pseudo 3D Voronoi cell of that node
    for (int i = 0; i < hexs.size(); ++i)
        // group hexs in vacuum
        if (hexs.get_marker(i) == TYPES.VACUUM) {
            for (int node : hexs[i])
                if (node >= node_min && node <= node_max) {
                    hexs.set_marker(i, 1 + node);
                    break;
                }
        // group hexs in bulk
        } else {
            for (int node : hexs[i])
                if (node >= node_min && node <= node_max) {
                    hexs.set_marker(i, -1 - node);
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

bool TetgenMesh::generate_hexahedra() {
    tethex::Mesh hexmesh;
    hexmesh.read_femocs(this);
    hexmesh.convert();
    hexmesh.export_femocs(this);

    group_hexahedra();
    calc_quad2hex2quad_mapping();
    nodes.calc_statistics();

    return 0;
}

int TetgenMesh::generate_surface(const string& cmd1, const string& cmd2) {
    TetgenMesh vacuum;
    vector<bool> tet_mask = vector_equal(tets.get_markers(), TYPES.VACUUM);

    // Transfer vacuum nodes and tetrahedra
    vacuum.nodes.copy(this->nodes);
    vacuum.tets.copy(this->tets, tet_mask);

    // calculate surface triangles
    int error_code = vacuum.recalc(cmd1);
    if (error_code) return error_code;

    // copy the triangles that are not on the simubox perimeter to input tetIO
    vacuum.tris.calc_appendices();
    const int n_surf_faces = tris.copy_surface(vacuum.tris, nodes.stat);

    // transfer elements and nodes to input
    nodes.transfer(false);
    tets.transfer(false);

    return calc_tet2tri_mapping(cmd2, n_surf_faces);
}

void TetgenMesh::generate_manual_surface() {
    const unsigned int n_elems = tets.size();
    const int max_surf_indx = nodes.indxs.surf_end;

    // booleans showing whether element i has exactly one face on the surface or not
    vector<bool> elem_on_surface; elem_on_surface.reserve(n_elems);
    // booleans showing whether node i is on the surface or not
    vector<bool> surf_locs(4);

    // Mark the elements that have exactly one face on the surface
    for (SimpleElement elem : tets) {
        for (unsigned int i = 0; i < 4; ++i)
            surf_locs[i] = (int)elem[i] <= max_surf_indx;
        elem_on_surface.push_back(vector_sum(surf_locs) == 3);
    }

//    // Reserve memory for surface faces
//    tris.init( vector_sum(elem_on_surface) + 2 );
//
//    // Make two big faces that pass the xy-min-max corners of surface
//    tris.append(SimpleFace(0, 1, 2));
//    tris.append(SimpleFace(0, 2, 3));

    // Reserve memory for surface faces
    tris.init( vector_sum(elem_on_surface) );

    // Generate the faces that separate material and vacuum
    // The faces are taken from the elements that have exactly one face on the surface
    for (unsigned int el = 0; el < n_elems; ++el)
        if (elem_on_surface[el]) {
            SimpleElement elem = tets[el];

            // Find the indices of nodes that are on the surface
            for (int i = 0; i < 4; ++i)
                surf_locs[i] = (int)elem[i] <= max_surf_indx;

            /* The possible combinations of surf_locs and n0,n1,n2:
             * surf_locs: 1110   1101   1011   0111
             *        n0: elem0  elem0  elem0  elem1
             *        n1: elem1  elem1  elem2  elem2
             *        n2: elem2  elem3  elem3  elem3   */
            int n0 = surf_locs[0] * elem[0] + (!surf_locs[0]) * elem[1];
            int n1 = (surf_locs[0] & surf_locs[1]) * elem[1] + (surf_locs[2] & surf_locs[3]) * elem[2];
            int n2 = (!surf_locs[3]) * elem[2] + surf_locs[3] * elem[3];
            tris.append(SimpleFace(n0, n1, n2));
        }
}

int TetgenMesh::tri2tet(const int tri, const int region) const {
    if (region == TYPES.VACUUM) {
        for (int tet : tris.to_tets(tri))
            if (tet >= 0 && tets.get_marker(tet) == TYPES.VACUUM)
                return tet;
    } else if (region == TYPES.BULK) {
        for (int tet : tris.to_tets(tri))
            if (tet >= 0 && tets.get_marker(tet) != TYPES.VACUUM)
                return tet;
    } else
        require(false, "Unimplemented region: " + d2s(region));

    return -1;
}

int TetgenMesh::quad2hex(const int quad, const int region) const {
    if (region == TYPES.VACUUM) {
        for (int hex : quads.to_hexs(quad))
            if (hex >= 0 && tets.get_marker(hexs.to_tet(hex)) == TYPES.VACUUM)
                return hex;
    } else if (region == TYPES.BULK) {
        for (int hex : quads.to_hexs(quad))
            if (hex >= 0 && tets.get_marker(hexs.to_tet(hex)) != TYPES.VACUUM)
                return hex;
    } else
        require(false, "Unimplemented region: " + d2s(region));

    return -1;
}

int TetgenMesh::calc_tet2tri_mapping(const string &cmd, int n_surf_faces) {
    tetIOout.deinitialize();
    tetIOout.initialize();

    // calculate the tetrahedron-triangle connectivity
    // the simubox boundary faces must also be calculated, no way to opt-out
    int error_code = recalc(cmd);
    if (error_code) return error_code;

    // the faces on the simubox sides are appended after the surface faces.
    // such property allows to remove the faces on the sides without affecting the tri2tet mapping.
    // such cleaning is useful to make other workflow faster.
    tris.init(n_surf_faces);
    for (int i = 0; i < n_surf_faces; ++i)
        tris.append(tris[i]);
    tris.transfer();

    vector<vector<int>> tet2tri_map(tets.size());
    for (int i = 0; i < n_surf_faces; ++i) {
        for (int tet : tris.to_tets(i))
            if (tet >= 0)
                tet2tri_map[tet].push_back(i);
    }
    tets.store_map(tet2tri_map);

    return 0;
}

void TetgenMesh::calc_quad2hex2quad_mapping() {
    const int n_quads = quads.size();
    const int n_hexs = hexs.size();

    vector<array<int,2>> quad2hex_map = vector<array<int,2>>(n_quads, {-1,-1});
    vector<vector<int>> hex2quad_map = vector<vector<int>>(n_hexs);

    for (int quad = 0; quad < n_quads; ++quad) {
        SimpleQuad squad = quads[quad];

        // loop through the tetrahedra that are connected to the quadrangle
        int region = 0;
        for (int tet : tris.to_tets(quads.to_tri(quad))) {
            if (tet < 0) continue;

            // loop through all the hexahedra connected to the tetrahedron
            for (int hex : tets.to_hexs(tet)) {

                // count for the # common nodes between quadrangle and hexahedron
                int n_common_nodes = 0;
                for (unsigned int node : hexs[hex])
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
    hexs.store_map(hex2quad_map);
}

int TetgenMesh::separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd) {
    vector<bool> tet_mask = vector_equal(tets.get_markers(), TYPES.VACUUM);
    vector<bool> hex_mask = vector_equal(hexs.get_markers(), TYPES.VACUUM);

    // Transfer vacuum nodes, tetrahedra, hexahedra and their markers
    vacuum.nodes.copy(this->nodes);
    vacuum.nodes.copy_markers(this->nodes);
    vacuum.tris.copy(this->tris);
    vacuum.tris.copy_markers(this->tris);
    vacuum.quads.copy(this->quads);
    vacuum.quads.copy_markers(this->quads);
    vacuum.tets.copy(this->tets, tet_mask);
    vacuum.tets.copy_markers(this->tets, tet_mask);
    vacuum.hexs.copy(this->hexs, hex_mask);
    vacuum.hexs.copy_markers(this->hexs, hex_mask);

    tet_mask.flip();
    hex_mask.flip();

    // Transfer bulk nodes, tetrahedra, hexahedra and their markers
    bulk.nodes.copy(this->nodes);
    bulk.nodes.copy_markers(this->nodes);
    bulk.tris.copy(this->tris);
    bulk.tris.copy_markers(this->tris);
    bulk.quads.copy(this->quads);
    bulk.quads.copy_markers(this->quads);
    bulk.tets.copy(this->tets, tet_mask);
    bulk.tets.copy_markers(this->tets, tet_mask);
    bulk.hexs.copy(this->hexs, hex_mask);
    bulk.hexs.copy_markers(this->hexs, hex_mask);

    return vacuum.recalc(cmd) + bulk.recalc(cmd);
}

void TetgenMesh::write_vtk(ofstream& out) {
    // k - write vtk, Q - quiet, I - suppresses iteration numbers,
    // F - suppress output of .face and .edge, E - suppress output of .ele
    const string cmd = "kIFEQ";
    const string path = "out/tetgenmesh";
    tetgenbehavior tetgenbeh;

    try {
        tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
        for (unsigned i = 0; i < path.size(); ++i)
            tetgenbeh.outfilename[i] = path[i];

        tetrahedralize(&tetgenbeh, &tetIOout, NULL);
    } catch (int e) {
        write_verbose_msg("TetgenMesh::write_vtk ended with error " + d2s(e));
    }
}

void TetgenMesh::write_bin(ofstream &out) const {
    const int n_nodes = nodes.stat.n_tetnode;
    const int n_tris = tris.size();
    const int n_tets = tets.size();

    // write Gmsh binary file header
    out << "$MeshFormat\n2.2 1 8\n";
    int one = 1;
    out.write((char*)&one, sizeof (int));
    out << "\n$EndMeshFormat\n$Nodes\n" << n_nodes << "\n";

    // write node coordinates
    for (int i = 0; i < n_nodes; ++i) {
        Point3 node = nodes[i];
        out.write ((char*)&i, sizeof (int));
        out.write ((char*)&node, sizeof (Point3));
    }

    out << "$\nEndNodes\n$Elements\n" << (n_tris + n_tets) << "\n";

    int header[3];
    int buffer[2] = {1, 0};

    // write header for triangles
    header[0]=Gmsh::triangle; header[1]=n_tris; header[2]=1;
    out.write ((char*)&header, sizeof (header));

    // write triangles
    for (int i = 0; i < n_tris; ++i, ++buffer[0]) {
        SimpleFace tri = tris[i];
        buffer[1] = tris.get_marker(i);
        out.write ((char*)&buffer, sizeof (buffer));
        out.write ((char*)&tri, sizeof (SimpleFace));
    }

    // write header for tetrahedra
    header[0]=Gmsh::tetrahedron; header[1]=n_tets; header[2]=1;
    out.write ((char*)&header, sizeof (header));

    // write tetrahedra
    for (int i = 0; i < n_tets; ++i, ++buffer[0]) {
        SimpleElement tet = tets[i];
        buffer[1] = tets.get_marker(i);
        out.write ((char*)&buffer, sizeof (buffer));
        out.write ((char*)&tet, sizeof (SimpleElement));
    }

    out << "\n$EndElements\n";
    out.close();
}

void TetgenMesh::write_msh(ofstream &out) const {
    // write Gmsh header
    FileWriter::write_msh(out);

    const int n_nodes = nodes.stat.n_tetnode;
    const int n_tris = tris.size();
    const int n_tets = tets.size();

    // write nodes
    out << "$Nodes\n" << n_nodes << "\n";

    for (size_t ver = 0; ver < n_nodes; ++ver)
        out << ver << " " << nodes[ver] << "\n";

    out << "$EndNodes\n$Elements\n" << (n_tris + n_tets) << "\n";

    int serial_nr = 1;

    // write triangles
    for (size_t i = 0; i < n_tris; ++i, ++serial_nr) {
        out << serial_nr << " "        // serial number of element
        << Gmsh::triangle << " 1 "     // Gmsh type of element & number of tags
        << tris.get_marker(i) << " "   // physical domain
        << tris[i] << "\n";
    }

    // write tetrahedra
    for (size_t i = 0; i < n_tets; ++i, ++serial_nr) {
        out << serial_nr << " "        // serial number of element
        << Gmsh::tetrahedron << " 1 "  // Gmsh type of element & number of tags
        << tets.get_marker(i) << " "   // physical domain
        << tets[i] << "\n";
    }

    out << "$EndElements\n";
}

void TetgenMesh::write_separate(const string& file_name, const int type) {
    if (not_write_time()) return;

    vector<bool> hex_mask;
    if (type == TYPES.VACUUM)
        hex_mask = vector_greater(hexs.get_markers(), 0);
    else
        hex_mask = vector_less(hexs.get_markers(), 0);

    TetgenMesh tempmesh;
    tempmesh.nodes.copy(this->nodes);
    tempmesh.nodes.copy_markers(this->nodes);
    tempmesh.nodes.transfer();
    tempmesh.hexs.copy(this->hexs, hex_mask);
    tempmesh.hexs.copy_markers(this->hexs, hex_mask);

    tempmesh.hexs.write(file_name);
}

void TetgenMesh::calc_pseudo_3D_vorocells(vector<vector<unsigned>>& cells, const bool vacuum) const {
    cells = vector<vector<unsigned>>(nodes.stat.n_tetnode);

    // the sign of hexahedron marker shows its location (>0 == vacuum, <0 == bulk)
    // and the (magnitude - 1) the index of tetrahedral node it is connected to
    int multiplier = 1;
    if (!vacuum) multiplier = -1;

    // find the pseudo Voronoi cell nodes for the tetrahedral nodes
    for (int hex = 0; hex < hexs.size(); ++hex) {
        int tetnode = multiplier * hexs.get_marker(hex);
        if (tetnode <= 0) continue;
        tetnode -= 1;

        for (int node : hexs[hex])
            if ( node != tetnode && nodes.get_marker(node) >= TYPES.EDGECENTROID )
                cells[tetnode].push_back(node);
    }
}

bool TetgenMesh::calc_ranks(vector<int>& ranks, const vector<vector<unsigned>>& nborlist) {
    const int n_nbor_layers = 4;  // number of nearest tetrahedra whose nodes will act as a seed
    const int n_nodes = nodes.size();
    const double max_rank = 100.0;
    const double eps = 0.01 * tets.stat.edgemin;

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

bool TetgenMesh::rank_and_mark_nodes() {
    const int n_nodes = nodes.size();
    const int min_rank = 30;
    int node;

    // Calculate neighbour list for nodes
    vector<vector<unsigned>> nborlist;
    tets.calc_nborlist(nborlist);

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

void TetgenMesh::mark_tets() {
    // Reserve memory for markers
    tets.init_markers(tets.size());

    // Locate all the elements
    for (SimpleElement elem : tets)
        tets.append_marker(locate_element(elem));
}

void TetgenMesh::mark_tris() {
    const double eps = 0.1 * tets.stat.edgemin;
    const int n_faces = tris.size();

    tris.init_markers(n_faces);
    nodes.calc_statistics();

    for (int i = 0; i < n_faces; ++i) {
        Point3 centre = tris.get_centroid(i);

        if (on_boundary(centre.x, nodes.stat.xmin, eps))
            tris.append_marker(TYPES.XMIN);
        else if (on_boundary(centre.x, nodes.stat.xmax, eps))
            tris.append_marker(TYPES.XMAX);
        else if (on_boundary(centre.y, nodes.stat.ymin, eps))
            tris.append_marker(TYPES.YMIN);
        else if (on_boundary(centre.y, nodes.stat.ymax, eps))
            tris.append_marker(TYPES.YMAX);
        else if (on_boundary(centre.z, nodes.stat.zmin, eps))
            tris.append_marker(TYPES.ZMIN);
        else if (on_boundary(centre.z, nodes.stat.zmax, eps))
            tris.append_marker(TYPES.ZMAX);
        else
            tris.append_marker(TYPES.SURFACE);
    }
}

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

bool TetgenMesh::mark_mesh() {
    if (rank_and_mark_nodes())
        return 1;

    mark_tets();
    return 0;
}

} /* namespace femocs */
