/*
 * TetgenCells.cpp
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#include <TetgenCells.h>
#include <float.h>

using namespace std;
namespace femocs {

// Initialize edge appending
const void TetgenEdges::init(const int N) {
    TetgenCells::init(N);
    writes->edgelist = new int[DIM * N];
}

// Append edge to mesh
const void TetgenEdges::append(const SimpleEdge &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->edgelist[i++] = node;
    i_cells++;
}

// Get i-th edge from the mesh
const SimpleCell<2> TetgenEdges::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleEdge(reads->edgelist[I], reads->edgelist[I+1]);
}

// Initialize face appending
const void TetgenFaces::init(const int N) {
    TetgenCells::init(N);
    writes->trifacelist = new int[DIM * N];
}

// Append face to mesh
const void TetgenFaces::append(const SimpleFace &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->trifacelist[i++] = node;
    i_cells++;
}

// Get i-th face from the mesh
const SimpleCell<3> TetgenFaces::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleFace(reads->trifacelist[I], reads->trifacelist[I+1], reads->trifacelist[I+2]);
}

// Initialize element appending
const void TetgenElements::init(const int N) {
    TetgenCells::init(N);
    writes->tetrahedronlist = new int[DIM * N];
}

// Append element to mesh
const void TetgenElements::append(const SimpleElement &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->tetrahedronlist[i++] = node;
    i_cells++;
}

// Get i-th element from the mesh
const SimpleCell<4> TetgenElements::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleElement(reads->tetrahedronlist[I], reads->tetrahedronlist[I+1],
            reads->tetrahedronlist[I+2], reads->tetrahedronlist[I+3]);
}

// Get indices of neighbouring elements of i-th element
const vector<int> TetgenElements::get_neighbours(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    require(reads->neighborlist, "Query from empty neighbour list!");

    const int I = DIM * i;
    const int* nborlist = reads->neighborlist;
    return vector<int> {nborlist[I+0], nborlist[I+1], nborlist[I+2], nborlist[I+3]};
}

// Initialize node appending
const void TetgenNodes::init(const int N) {
    TetgenCells::init(N);
    writes->pointlist = new double[n_coordinates * N];
}

// Append node to mesh
const void TetgenNodes::append(const Point3 &point) {
    require(i_cells < *n_cells_w, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_cells;
    for (double node : point)
        writes->pointlist[i++] = node;
    i_cells++;
}

// Get i-th node from the mesh
const SimpleCell<1> TetgenNodes::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return SimpleNode(i);
}

// Modify the coordinates of i-th node
const void TetgenNodes::set_node(const int i, const Point3 &point) {
    require(i >= 0 && i < get_n_nodes(), "Index out of bounds: " + to_string(i));
    int I = n_coordinates * i;
    for (double node : point)
        reads->pointlist[I++] = node;
}

// Return the coordinates of i-th node as a 3D vector
const Vec3 TetgenNodes::get_vec(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int n = n_coordinates * i;
    return Vec3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
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
    size_t n_nodes = size();
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

// Copy the nodes from another mesh
const void TetgenNodes::copy(const TetgenNodes& nodes, const vector<bool>& mask) {
    const int n_nodes = nodes.size();

    // In case of empty or non-aligned mask, copy all the nodes
    if (n_nodes != mask.size()) {
        init(n_nodes);
        for (int i = 0; i < n_nodes; ++i)
            append(nodes[i]);
        i_cells = n_nodes;

    // In case of aligned mask, copy only the nodes specified by the mask
    } else {
        const int n_mask = vector_sum(mask);
        init(n_mask);
        for (int i = 0; i < n_nodes; ++i)
            if (mask[i])
                append(nodes[i]);
        i_cells = n_mask;
    }
}

// Write node data to file
const void TetgenNodes::write(const string &file_name) const {
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
const void TetgenNodes::write_xyz(const string &file_name) const {
    const int n_nodes = size();
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
