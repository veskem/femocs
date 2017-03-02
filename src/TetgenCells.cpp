/*
 * TetgenCells.cpp
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#include "TetgenCells.h"
#include <float.h>

using namespace std;
namespace femocs {

/* =====================================================================
 *  =========================== TetgenNodes ===========================
 * ===================================================================== */

// Copy the nodes from write to read buffer
void TetgenNodes::recalc() {
    TetgenCells::recalc();
    reads->pointlist = new double[n_coordinates * i_cells];
    for (int i = 0; i < n_coordinates * i_cells; ++i)
        reads->pointlist[i] = writes->pointlist[i];
}

// Initialize node appending
void TetgenNodes::init(const int N) {
    TetgenCells::init(N);
    writes->pointlist = new double[n_coordinates * N];
}

// Append node to mesh
void TetgenNodes::append(const Point3 &point) {
    require(i_cells < *n_cells_w, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_cells;
    for (double node : point)
        writes->pointlist[i++] = node;
    i_cells++;
}

// Get i-th node from the mesh
SimpleCell<1> TetgenNodes::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return SimpleNode(i);
}

// Return the coordinates of i-th node as a 3D vector
Vec3 TetgenNodes::get_vec(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int n = n_coordinates * i;
    return Vec3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
}

// Modify the coordinates of i-th node
void TetgenNodes::set_node(const int i, const Point3 &point) {
    require(i >= 0 && i < get_n_nodes(), "Index out of bounds: " + to_string(i));
    int I = n_coordinates * i;
    for (double node : point)
        reads->pointlist[I++] = node;
}

void TetgenNodes::save_indices(const int n_surf, const int n_bulk, const int n_vacuum) {
    indxs.surf_start = 0;
    indxs.surf_end = indxs.surf_start + n_surf - 1;
    indxs.bulk_start = indxs.surf_end + 1;
    indxs.bulk_end = indxs.bulk_start + n_bulk - 1;
    indxs.vacuum_start = indxs.bulk_end + 1;
    indxs.vacuum_end = indxs.vacuum_start + n_vacuum - 1;
    indxs.tetgen_start = indxs.vacuum_end + 1;
    indxs.tetgen_end = -1;
}

// Store the locations of different kinds of nodes that were produced while splitting tetrahedra into hexahedra
void TetgenNodes::save_hex_indices(const vector<int>& n_nodes) {
    require(n_nodes.size() == 4, "Invalid indices!");

    stat.n_tetnode = n_nodes[0];
    stat.n_midedge = n_nodes[1];
    stat.n_midface = n_nodes[2];
    stat.n_midtet  = n_nodes[3];

    indxs.tetnode_start = 0;
    indxs.tetnode_end   = indxs.tetnode_start + n_nodes[0] - 1;
    indxs.midedge_start = indxs.tetnode_end + 1;
    indxs.midedge_end   = indxs.midedge_start + n_nodes[1] - 1;
    indxs.midface_start = indxs.midedge_end + 1;
    indxs.midface_end   = indxs.midface_start + n_nodes[2] - 1;
    indxs.midtet_start  = indxs.midface_end + 1;
    indxs.midtet_end    = indxs.midtet_start  + n_nodes[3] - 1;
}

void TetgenNodes::init_statistics() {
    stat.n_bulk = stat.n_surface = stat.n_vacuum = 0;
    stat.xmin = stat.ymin = stat.zmin = DBL_MAX;
    stat.xmax = stat.ymax = stat.zmax =-DBL_MAX;
}

void TetgenNodes::calc_statistics() {
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
void TetgenNodes::copy(const TetgenNodes& nodes, const vector<bool>& mask) {
    const int n_nodes = nodes.size();
    copy_statistics(nodes);

    // In case of empty or non-aligned mask, copy all the nodes
    if (n_nodes != mask.size()) {
        init(n_nodes);
        for (int i = 0; i < n_nodes; ++i)
            append(nodes[i]);

    // In case of aligned mask, copy only the nodes specified by the mask
    } else {
        const int n_mask = vector_sum(mask);
        init(n_mask);
        for (int i = 0; i < n_nodes; ++i)
            if (mask[i]) append(nodes[i]);
    }
}

void TetgenNodes::copy_statistics(const TetgenNodes& n) {
    require(this->indxs.size() == n.indxs.size() , "Incompatible indices!");
    for (int i = 0; i < indxs.size(); ++i)
        (&this->indxs.surf_start)[i] = (&n.indxs.surf_start)[i];

    stat.n_tetnode = n.stat.n_tetnode;
    stat.n_midedge = n.stat.n_midedge;
    stat.n_midface = n.stat.n_midface;
    stat.n_midtet  = n.stat.n_midtet;
}

// Transform nodes into Deal.II format
vector<dealii::Point<3>> TetgenNodes::export_dealii() const {
    vector<dealii::Point<3>> nodes; nodes.reserve(get_n_nodes());
    for (Point3 p : *this)
        nodes.push_back( dealii::Point<3>(p.x, p.y, p.z) );
    return nodes;
}

// Write node data to file
void TetgenNodes::write(const string &file_name) const {
    if (!MODES.WRITEFILE) return;

    string file_type = get_file_type(file_name);
    require(file_type == "xyz" || file_type == "vtk", "Unknown file type: " + file_type);

    const int celltype = 1;     // 1-vertex, 3-line, 5-triangle, 10-tetrahedron

    if (file_type == "xyz")
        write_xyz(file_name);
    else
        write_vtk(file_name, celltype);
}

// Write node data to .xyz file
void TetgenNodes::write_xyz(const string &file_name) const {
    const int n_nodes = size();
    const int n_markers = get_n_markers();

    expect(n_nodes > 0, "Zero nodes detected!");

    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << fixed;
    out_file << n_nodes << "\n";
    out_file << "Mesh nodes properties=id:R:1:pos:R:3:marker:R:1\n";

    if (n_nodes == n_markers)
        for (int i = 0; i < n_nodes; ++i)
            out_file << i << " " << get_node(i) << " " << get_marker(i) << endl;
    else
        for (int i = 0; i < n_nodes; ++i)
            out_file << i << " " << get_node(i) << " " << "-1" << endl;

    out_file.close();
}

/* =====================================================================
 *  =========================== TetgenEdges ===========================
 * ===================================================================== */

// Copy the nodes from write to read buffer
void TetgenEdges::recalc() {
    TetgenCells::recalc();
    reads->edgelist = new int[DIM * i_cells];
    for (int i = 0; i < DIM * i_cells; ++i)
        reads->edgelist[i] = writes->edgelist[i];
}

// Initialize edge appending
void TetgenEdges::init(const int N) {
    TetgenCells::init(N);
    writes->edgelist = new int[DIM * N];
}

// Append edge to mesh
void TetgenEdges::append(const SimpleCell<2> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->edgelist[i++] = node;
    i_cells++;
}

// Get i-th edge from the mesh
SimpleCell<2> TetgenEdges::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleEdge(reads->edgelist[I], reads->edgelist[I+1]);
}

/* =====================================================================
 *  =========================== TetgenFaces ===========================
 * ===================================================================== */

// Copy the nodes from write to read buffer
void TetgenFaces::recalc() {
    TetgenCells::recalc();
    reads->trifacelist = new int[DIM * i_cells];
    for (int i = 0; i < DIM * i_cells; ++i)
        reads->trifacelist[i] = writes->trifacelist[i];
}

// Initialize face appending
void TetgenFaces::init(const int N) {
    TetgenCells::init(N);
    writes->trifacelist = new int[DIM * N];
}

// Append face to mesh
void TetgenFaces::append(const SimpleCell<3> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->trifacelist[i++] = node;
    i_cells++;
}

// Get i-th face from the mesh
SimpleCell<3> TetgenFaces::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleFace(reads->trifacelist[I], reads->trifacelist[I+1], reads->trifacelist[I+2]);
}

// Delete the faces on the sides of simulation cell
//void TetgenFaces::clean_sides(const TetgenNodes::Stat& stat) {
void TetgenFaces::clean_sides(const Medium::Sizes& stat) {
    const double eps = 0.1;
    const int n_faces = size();
    vector<SimpleFace> faces; faces.reserve(n_faces);

    // Make a copy of faces not on the sides of simulation cell
    for (int i = 0; i < n_faces; ++i) {
        const Point3 centroid = get_centroid(i);
        const bool side_x = on_boundary(centroid.x, stat.xmin, stat.xmax, eps);
        const bool side_y = on_boundary(centroid.y, stat.ymin, stat.ymax, eps);
        const bool side_z = on_boundary(centroid.z, stat.zmin, stat.zmax, eps);
        if (!(side_x || side_y || side_z))
            faces.push_back(get_cell(i));
    }

    // Store the surface faces
    init(faces.size());
    for (SimpleFace face : faces)
        append(face);

    calc_norms();
}

// Calculate the norms and areas for all the triangles
void TetgenFaces::calc_norms() {
    const int n_faces = size();
    areas.clear(); areas.reserve(n_faces);
    norms.clear(); norms.reserve(n_faces);
    
    for (SimpleFace sface : *this) {
        Vec3 v0 = get_vec(sface[0]);
        Vec3 v1 = get_vec(sface[1]);
        Vec3 v2 = get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 n = e1.crossProduct(e2);
        
        areas.push_back(n.norm() * 0.5);
        norms.push_back(n.normalize());
    }
}

// Return the normal of i-th triangle
Vec3 TetgenFaces::get_norm(const int i) const {
    require(i >= 0 && i < norms.size(), "Invalid index: " + to_string(i));
    return norms[i];
}

// Return the area of i-th triangle
double TetgenFaces::get_area(const int i) const {
    require(i >= 0 && i < areas.size(), "Invalid index: " + to_string(i));
    return areas[i];
}
    
/* =====================================================================
 *  ========================== TetgenElements =========================
 * ===================================================================== */

// Copy the nodes from write to read buffer
void TetgenElements::recalc() {
    *n_cells_r = *n_cells_w;
    i_cells = *n_cells_w;
    reads->tetrahedronlist = new int[DIM * i_cells];
    for (int i = 0; i < DIM * i_cells; ++i)
        reads->tetrahedronlist[i] = writes->tetrahedronlist[i];
}

// Initialize element appending
void TetgenElements::init(const int N) {
    TetgenCells::init(N);
    writes->tetrahedronlist = new int[DIM * N];
}

// Append element to mesh
void TetgenElements::append(const SimpleCell<4> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->tetrahedronlist[i++] = node;
    i_cells++;
}

// Get i-th element from the mesh
SimpleCell<4> TetgenElements::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    const int I = DIM * i;
    return SimpleElement(reads->tetrahedronlist[I], reads->tetrahedronlist[I+1],
            reads->tetrahedronlist[I+2], reads->tetrahedronlist[I+3]);
}

// Get indices of neighbouring elements of i-th element
vector<int> TetgenElements::get_neighbours(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    require(reads->neighborlist, "Query from empty neighbour list!");

    const int I = DIM * i;
    const int* nborlist = reads->neighborlist;
    return vector<int> {nborlist[I+0], nborlist[I+1], nborlist[I+2], nborlist[I+3]};
}

/* =====================================================================
 *  ============================ Hexahedra ============================
 * ===================================================================== */

// Get number of hexahedra in mesh
int Hexahedra::size() const {
    return hexs.size();
}

// Initialize hexahedron appending
void Hexahedra::init(const int N) {
    TetgenCells::init(N);
    hexs.clear();
    hexs.reserve(N);
}

// Append hexahedron to mesh
void Hexahedra::append(const SimpleCell<8> &cell) {
    expect(hexs.size() < hexs.capacity(), "Allocated size of cells exceeded!");
    hexs.push_back(cell);
    i_cells++;
}

// Get i-th element from the mesh
SimpleCell<8> Hexahedra::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return hexs[i];
}

// Transform hexahedra into Deal.II format
vector<dealii::CellData<3>> Hexahedra::export_dealii() const {
    const int n_elems = size();
    std::vector<dealii::CellData<3>> elems; elems.reserve(n_elems);

    // loop through all the hexahedra
    for (SimpleHex hex : *this) {
        elems.push_back(dealii::CellData<3>());
        for (unsigned int v = 0; v < DIM; ++v)
            elems.back().vertices[v] = hex[v];
    }
    return elems;
}

/* =====================================================================
 *  ============================ Voronoi ============================
 * ===================================================================== */

Vec3 VoronoiFace::norm() {
    Vec3 edge1 = verts[1] - verts[0];
    Vec3 edge2 = verts[2] - verts[0];
    return edge1.crossProduct(edge2).normalize();
}

Vec3 VoronoiFace::area() {
    const int n_nodes = size();
    calc_nodes();

    Vec3 total(0.0);
    for (int i = 0; i < n_nodes; ++i) {
        int j = i+1;
        if (j == n_nodes) j = 0;
        total += verts[i].crossProduct(verts[j]);
    }

    Vec3 nrm = norm();
    double area = fabs(0.5 * total.dotProduct(nrm));

    return nrm * area;
}

Vec3 VoronoiFace::centroid() {
    calc_nodes();
    Vec3 centroid(0.0);

    for (Vec3& v : verts)
        centroid += v;
    centroid *= (1.0 / size());

    return centroid;
}

// Calculate the unique node that is associated with the edge
int VoronoiFace::get_node(const int edge) {
    static int node = -1;

    const int v1 = data->vedgelist[edge].v1;
    if (node != v1) node = v1;
    else node = data->vedgelist[edge].v2;

    return node;
}

// Transform the node data from tetgenio into easily accessible form
void VoronoiFace::calc_nodes() {
    if (verts.size() > 0) return;
    verts.reserve(size());

    for (int edge : *this) {
        int n = 3 * get_node(edge);
        verts.push_back( Vec3(data->vpointlist[n+0], data->vpointlist[n+1], data->vpointlist[n+2]) );
    }
}

// Return the neighbouring cell for the caller cell
int VoronoiFace::nborcell(const int caller_id) {
    const int c1 = data->vfacetlist[id].c1;
    if (c1 == caller_id)
        return data->vfacetlist[id].c2;
    return c1;
}

void VoronoiCells::write_cells(ofstream& out) const {
    // Get total number of faces and their vertices
    size_t n_faces = 0, n_verts = 0;
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell) {
            n_verts += face.size();
            n_faces++;
        }

    // Write the Voronoi cells as polygons
    out << "\nCELLS " << n_faces << " " << n_faces + n_verts << "\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << face.size() << " " << face << endl;

    // Output cell types
    out << "\nCELL_TYPES " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; ++i)
        out << celltype << "\n";

    // Output scalar data associated with Voronoi cells
    out << "\nCELL_DATA " << n_faces << "\n";

    // Write the IDs of cells
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << cell.id << "\n";

    // Write the markers of cells
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (VoronoiCell cell : *this)
        for (VoronoiFace face : cell)
            out << get_marker(cell.id) << "\n";
}

// Calculate the face neighbours - faces, that share the same edge
void VoronoiFaces::calc_neighbours() {
    const int n_faces = size();
    const int n_edges = tetio->numberofvedges;

    // Produce edge to face mapping
    vector<vector<int>> edge2face(n_edges);
    for (int i = 0; i < n_edges; ++i)
        edge2face[i].reserve(10);

    for (VoronoiFace face : *this)
        for (int i : face)
            edge2face[i].push_back(face.id);

    // Initialize face neighbour list
    neighbours.clear();
    neighbours.resize(n_faces);
    for (int i = 0; i < n_faces; ++i)
        neighbours[i].reserve(10);

    // Convert edge to face list to face neighbour list
    for (vector<int> faces : edge2face)
        for (int i : faces)
            for (int j : faces) {
                if (i == j) continue;
                neighbours[i].push_back(j);
            }
}

void VoronoiFaces::calc_centroids() {
    const int n_faces = size();
    centroids.clear();
    centroids.reserve(n_faces);

    for (VoronoiFace face : *this) {
        if (face.size() > 0)
            centroids.push_back(face.centroid());
        else
            centroids.push_back(Vec3());
    }
}

void VoronoiFaces::write_cells(ofstream& out) const {
    // Get total number of faces and their vertices
    size_t n_faces = 0, n_verts = 0;
    for (VoronoiFace face : *this) {
        n_verts += face.size();
        n_faces += face.size() > 0;
    }

    // Write the Voronoi faces as polygons
    out << "\nCELLS " << n_faces << " " << n_faces + n_verts << "\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << face.size() << " " << face << endl;

    // Output cell types
    out << "\nCELL_TYPES " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; ++i)
        out << celltype << "\n";

    // Output scalar data associated with Voronoi faces
    out << "\nCELL_DATA " << n_faces << "\n";

    // write IDs of faces
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << face.id << "\n";

    // write markers of faces
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (VoronoiFace face : *this)
        if (face.size() > 0)
            out << get_marker(face.id) << "\n";
}

} /* namespace femocs */
