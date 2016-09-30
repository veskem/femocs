/*
 * TetgenCells.cpp
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#include <TetgenCells.h>
#include <fstream>

using namespace std;
namespace femocs {

// Get i-th edge from the mesh
const SimpleEdge TetgenEdges::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleEdge(nodes_r[I], nodes_r[I+1]);
}

// Get i-th face from the mesh
const SimpleFace TetgenFaces::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleFace(nodes_r[I], nodes_r[I+1], nodes_r[I+2]);
}

// Get i-th element from the mesh
const SimpleElement TetgenElements::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleElement(nodes_r[I], nodes_r[I+1], nodes_r[I+2], nodes_r[I+3]);
}

// Get indices of neighbouring elements of i-th element
const vector<int> TetgenElements::get_neighbours(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    require(neighborlist, "Query from empty neighbour list!");

    const int I = DIM * i;
    return vector<int> {neighborlist[I+0], neighborlist[I+1], neighborlist[I+2], neighborlist[I+3]};
}

// Get i-th node from the mesh
const SimpleNode TetgenNodes::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    return SimpleNode(i);
}

// Return the coordinates of i-th node as a 3D vector
const Vec3 TetgenNodes::get_vec(const int i) const {
    require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
    const int n = n_coordinates * i;
    return Vec3(points_r[n+0], points_r[n+1], points_r[n+2]);
}

// Initialize node appending
const void TetgenNodes::init_nodes(const int N) {
    require(N > 0, "Invalid number of nodes: " + to_string(N));
    i_cells = 0;
    *n_points_w = N;

    points_w = new double[n_coordinates * N];
}

// Append node to mesh
const void TetgenNodes::add_node(const Point3 &point) {
    require(i_cells < *n_points_w, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_cells;
    for (double node : point)
        points_w[i++] = node;
    i_cells++;
}

// Write node data to file
const void TetgenNodes::write(const string &file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "xyz" || file_type == "vtk", "Unknown file type: " + file_type);

    const int celltype = 1;     // 1-vertex, 3-line, 5-triangle, 10-tetrahedron

    if (file_type == "xyz")
        write_xyz(file_name);
    else
        write_vtk(file_name, celltype);
}

// Write node data to .xyz file
const void TetgenNodes::write_xyz(const string &file_name) {
    const int n_nodes = get_n_nodes();
    const int n_markers = get_n_markers();

    expect(n_nodes > 0, "Zero nodes detected!");

    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << n_nodes << "\n";
    out_file << "Mesh nodes: id x y z marker\n";

    if (n_nodes == n_markers)
        for (int i = 0; i < n_nodes; ++i)
            out_file << i << " " << get_node(i) << " " << get_marker(i) << endl;
    else
        for (int i = 0; i < n_nodes; ++i)
            out_file << i << " " << get_node(i) << " " << "-1" << endl;

    out_file.close();
}

} /* namespace femocs */
