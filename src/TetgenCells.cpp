/*
 * TetgenCells.cpp
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#include <TetgenCells.h>
#include <fstream>
#include <float.h>

using namespace std;
namespace femocs {

template<size_t dim>
const void TetgenCells<dim>::write_vtk(const string &file_name, const int celltype) {
    const int n_markers = get_n_markers();
    const int n_nodes = get_n_nodes();
    const int n_cells = *n_cells_r;

    expect(n_nodes > 0, "Zero nodes detected!");

    std::ofstream out(file_name.c_str());
    require(out, "Can't open a file " + file_name);

    out.setf(std::ios::scientific);
    out.precision(8);

    out << "# vtk DataFile Version 3.0\n";
    out << "# TetgenCells data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the nodes
    out << "POINTS " << n_nodes << " double\n";
    for (size_t ver = 0; ver < n_nodes; ++ver)
        out << get_node(ver) << "\n";

    // Output the cells (tetrahedra, triangles, edges or vertices)
    out << "\nCELLS " << n_cells << " " << n_cells * (dim + 1) << "\n";
    for (size_t el = 0; el < n_cells; ++el)
        out << dim << " " << get_cell(el) << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (size_t el = 0; el < n_cells; ++el)
        out << celltype << "\n";

    // Output cell markers
    if ((n_markers > 0) && (n_markers == n_cells)) {
        out << "\nCELL_DATA " << n_cells << "\n";
        out << "SCALARS Cell_markers int\nLOOKUP_TABLE default\n";
        for (size_t el = 0; el < n_cells; ++el)
            out << get_marker(el) << "\n";
    }
}

// Get i-th edge from the mesh
const SimpleCell<2> TetgenEdges::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleEdge(nodes_r[I], nodes_r[I+1]);
}

// Get i-th face from the mesh
const SimpleCell<3> TetgenFaces::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleFace(nodes_r[I], nodes_r[I+1], nodes_r[I+2]);
}

// Get i-th element from the mesh
const SimpleCell<4> TetgenElements::get_cell(const int i) const {
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
const SimpleCell<1> TetgenNodes::get_cell(const int i) const {
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
const void TetgenNodes::init(const int N) {
    require(N > 0, "Invalid number of nodes: " + to_string(N));
    i_cells = 0;
    *n_points_w = N;
    points_w = new double[n_coordinates * N];
}

// Append node to mesh
const void TetgenNodes::append(const Point3 &point) {
    require(i_cells < *n_points_w, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_cells;
    for (double node : point)
        points_w[i++] = node;
    i_cells++;
}

const void TetgenNodes::save_indices(const int n_surf, const int n_bulk, const int n_vacuum) {
    indxs.surf_start = 0;
    indxs.surf_end = indxs.surf_start + n_surf - 1;
    indxs.bulk_start = indxs.surf_end + 1;
    indxs.bulk_end = indxs.bulk_start + n_bulk - 1;
    indxs.vacuum_start = indxs.bulk_end + 1;
    indxs.vacuum_end = indxs.vacuum_start + n_vacuum - 1;
    indxs.tetgen_start = indxs.vacuum_end + 1;
}

const void TetgenNodes::init_statistics() {
    stat.n_bulk = stat.n_surface = stat.n_vacuum = 0;
    stat.xmin = stat.ymin = stat.zmin = DBL_MAX;
    stat.xmax = stat.ymax = stat.zmax =-DBL_MAX;
}

const void TetgenNodes::calc_statistics() {
    init_statistics();
    size_t n_nodes = get_n_nodes();
    size_t n_markers = get_n_markers();

    // Find the min and max coordinates of all nodes
    for (int i = 0; i < n_nodes; ++i) {
        Point3 point = get_node(i);
        stat.xmax = max(stat.xmax, point.x);
        stat.xmin = min(stat.xmin, point.x);
        stat.ymax = max(stat.ymax, point.y);
        stat.ymin = min(stat.ymin, point.y);
        stat.zmax = max(stat.zmax, point.z);
        stat.zmin = min(stat.zmin, point.z);
    }

    // Find the number of nodes in various regions
    if (n_markers > 0)
        // Loop through all the node markers
        for (int marker : markers) {
            if (marker == TYPES.BULK)
                stat.n_bulk++;
            else if (marker == TYPES.VACUUM)
                stat.n_vacuum++;
            else if (marker == TYPES.SURFACE)
                stat.n_surface++;
        }
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
