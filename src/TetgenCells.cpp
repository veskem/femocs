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

void TetgenNodes::transfer(const bool write2read) {
    if (reads == writes) return;
    TetgenCells::transfer(write2read);
    init_statistics();
    if (write2read) {
        if (reads->pointlist != (REAL*) NULL)
            delete[] reads->pointlist;
        reads->pointlist = new REAL[n_coordinates * i_cells];
        for (int i = 0; i < n_coordinates * i_cells; ++i)
            reads->pointlist[i] = writes->pointlist[i];
    } else {
        if (writes->pointlist != (REAL*) NULL)
            delete[] writes->pointlist;
        writes->pointlist = new REAL[n_coordinates * i_cells];

        if (writes->pointmarkerlist != (int*) NULL)
            delete[] writes->pointmarkerlist;
        writes->pointmarkerlist = new int[i_cells];

        for (int i = 0; i < n_coordinates * i_cells; ++i)
            writes->pointlist[i] = reads->pointlist[i];
    }
}

void TetgenNodes::init(const int N) {
    TetgenCells::init(N);
    init_statistics();
    if (writes->pointlist != (REAL*) NULL)
        delete[] writes->pointlist;
    writes->pointlist = new REAL[n_coordinates * N];
}

void TetgenNodes::append(const Point3 &point) {
    require(i_cells < *n_cells_w, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_cells;
    for (double node : point)
        writes->pointlist[i++] = node;
    i_cells++;
}

SimpleCell<1> TetgenNodes::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    return SimpleNode(i);
}

Vec3 TetgenNodes::get_vec(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    const int n = n_coordinates * i;
    return Vec3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
}

void TetgenNodes::set_node(const int i, const Point3 &point) {
    require(i >= 0 && i < get_n_nodes(), "Index out of bounds: " + d2s(i));
    int I = n_coordinates * i;
    for (double node : point)
        reads->pointlist[I++] = node;
}

void TetgenNodes::set_boundary(const int i, const int value) {
    require(i >= 0 && i < get_n_nodes(), "Index out of bounds: " + d2s(i));
    writes->pointmarkerlist[i] = value;
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

void TetgenNodes::save_hex_indices(const vector<int>& n_nodes) {
    require(n_nodes.size() == 4, "Invalid size of indices: " + d2s(n_nodes.size()));

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
    stat.xbox = stat.ybox = stat.zbox = 0;
}

void TetgenNodes::calc_statistics() {
    init_statistics();
    int n_nodes = size();
    int n_markers = get_n_markers();

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
    stat.xbox = stat.xmax - stat.xmin;
    stat.ybox = stat.ymax - stat.ymin;
    stat.zbox = stat.zmax - stat.zmin;

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

void TetgenNodes::copy(const TetgenNodes& nodes, const vector<bool>& mask) {
    const size_t n_nodes = nodes.size();
    copy_statistics(nodes);

    // In case of empty or non-aligned mask, copy all the nodes
    if (n_nodes != mask.size()) {
        init(n_nodes);
        for (size_t i = 0; i < n_nodes; ++i)
            append(nodes[i]);

    // In case of aligned mask, copy only the nodes specified by the mask
    } else {
        const int n_mask = vector_sum(mask);
        init(n_mask);
        for (size_t i = 0; i < n_nodes; ++i)
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

vector<dealii::Point<3>> TetgenNodes::export_dealii() const {
    vector<dealii::Point<3>> nodes; nodes.reserve(get_n_nodes());
    for (Point3 p : *this)
        nodes.push_back( dealii::Point<3>(p.x, p.y, p.z) );
    return nodes;
}

void TetgenNodes::write_xyz(ofstream &out) const {
    const int n_nodes = size();
    const int n_markers = get_n_markers();

    FileWriter::write_xyz(out); // write start of the header
    out << "properties=id:R:1:pos:R:3:marker:R:1\n";

    if (n_nodes == n_markers)
        for (int i = 0; i < n_nodes; ++i)
            out << i << " " << get_node(i) << " " << get_marker(i) << endl;
    else
        for (int i = 0; i < n_nodes; ++i)
            out << i << " " << get_node(i) << " " << "-1" << endl;
}

/* =====================================================================
 *  =========================== TetgenEdges ===========================
 * ===================================================================== */

void TetgenEdges::transfer(const bool write2read) {
    if (reads == writes) return;
    TetgenCells::transfer(write2read);
    if (write2read) {
        if (reads->edgelist != (int *) NULL)
            delete[] reads->edgelist;
        reads->edgelist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            reads->edgelist[i] = writes->edgelist[i];
    } else {
        if (writes->edgelist != (int *) NULL)
            delete[] writes->edgelist;
        writes->edgelist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            writes->edgelist[i] = reads->edgelist[i];
    }
}

void TetgenEdges::init(const int N) {
    TetgenCells::init(N);
    if (writes->edgelist != (int *) NULL)
        delete[] writes->edgelist;
    writes->edgelist = new int[DIM * N];
}

void TetgenEdges::append(const SimpleCell<2> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->edgelist[i++] = node;
    i_cells++;
}

SimpleCell<2> TetgenEdges::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    const int I = DIM * i;
    return SimpleEdge(reads->edgelist[I], reads->edgelist[I+1]);
}

void TetgenEdges::init_statistics() {
    stat.edgemin =  DBL_MAX;
    stat.edgemax = -DBL_MAX;
}

void TetgenEdges::calc_statistics() {
    init_statistics();

    // loop through all the triangles
    for (SimpleEdge edge : *this) {
        // read the coordinates of triangle nodes
        Point3 n1 = get_node(edge[0]);
        Point3 n2 = get_node(edge[1]);

        // calculate squared lengts of the edge
        const double e = n1.distance2(n2);

        // store min and max length
        stat.edgemin = min(stat.edgemin, e);
        stat.edgemax = max(stat.edgemax, e);
    }

    // squared length -> length
    stat.edgemin = sqrt(stat.edgemin);
    stat.edgemax = sqrt(stat.edgemax);
}

void TetgenEdges::clean_sides(const Medium::Sizes& stat) {
    const double eps = 0.01 * this->stat.edgemin;
    const int n_faces = size();
    vector<SimpleEdge> edges; edges.reserve(n_faces);

    // Make a copy of edges on the sides of simulation cell
    for (int i = 0; i < n_faces; ++i) {
        const Point3 centroid = get_centroid(i);
        const bool side_x = on_boundary(centroid.x, stat.xmin, stat.xmax, eps);
        const bool side_y = on_boundary(centroid.y, stat.ymin, stat.ymax, eps);
        const bool side_z = on_boundary(centroid.z, stat.zmin, stat.zmax, eps);
        if (!side_z && (side_x || side_y) && !(side_x && side_y))
            edges.push_back(get_cell(i));
    }

    // Store the edges that are on the simulation cell boundary
    init(edges.size());
    for (SimpleEdge edge : edges)
        append(edge);
}

/* =====================================================================
 *  =========================== TetgenFaces ===========================
 * ===================================================================== */

void TetgenFaces::transfer(const bool write2read) {
    if (reads == writes) return;
    TetgenCells::transfer(write2read);
    if (write2read) {
        if (reads->trifacelist != (int *) NULL)
            delete[] reads->trifacelist;
        reads->trifacelist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            reads->trifacelist[i] = writes->trifacelist[i];
    } else {
        if (writes->trifacelist != (int *) NULL)
            delete[] writes->trifacelist;
        writes->trifacelist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            writes->trifacelist[i] = reads->trifacelist[i];
    }
}

void TetgenFaces::init(const int N) {
    TetgenCells::init(N);
    if (writes->trifacelist != (int *) NULL)
        delete[] writes->trifacelist;
    writes->trifacelist = new int[DIM * N];
}

void TetgenFaces::init_with_map(const int N) {
    require(N >= 0, "Invalid number of cells: " + d2s(N));
    i_cells = 0;
    *n_cells_r = N;

    if (reads->trifacelist != (int *) NULL)
        delete[] reads->trifacelist;
    reads->trifacelist = new int[DIM * N];

    if (reads->face2tetlist != (int *) NULL)
        delete[] reads->face2tetlist;
    reads->face2tetlist = new int[n_tets_per_tri * N];
}

void TetgenFaces::append(const SimpleCell<3> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->trifacelist[i++] = node;
    i_cells++;
}

void TetgenFaces::append(const SimpleCell<3> &cell, int tet1, int tet2) {
    require(i_cells < *n_cells_r, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        reads->trifacelist[i++] = node;

    i = n_tets_per_tri * i_cells;
    reads->face2tetlist[i++] = tet1;
    reads->face2tetlist[i] = tet2;

    i_cells++;
}

SimpleCell<3> TetgenFaces::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    const int I = DIM * i;
    return SimpleFace(reads->trifacelist[I], reads->trifacelist[I+1], reads->trifacelist[I+2]);
}

int TetgenFaces::copy_surface(const TetgenFaces& faces, const TetgenNodes::Stat& stat) {
    const double eps = 0.01 * faces.stat.edgemin;
    int n_faces = faces.size();
    vector<bool> tri_mask(n_faces);

    // Make a copy of faces not on the sides of simulation cell
    for (int i = 0; i < n_faces; ++i) {
        const Point3 centroid = faces.get_centroid(i);
        const bool side_x = on_boundary(centroid.x, stat.xmin, stat.xmax, eps);
        const bool side_y = on_boundary(centroid.y, stat.ymin, stat.ymax, eps);
        const bool side_z = on_boundary(centroid.z, stat.zmin, stat.zmax, eps);
        tri_mask[i] = !(side_x || side_y || side_z);
    }

    // Store the surface faces
    int n_surface_faces = vector_sum(tri_mask);
    init(n_surface_faces);
    for (int i = 0; i < n_faces; ++i)
        if (tri_mask[i])
            append(faces[i]);

    calc_appendices();
    return n_surface_faces;
}

void TetgenFaces::calc_appendices() {
    const int n_faces = size();
    areas.clear(); areas.reserve(n_faces);
    norms.clear(); norms.reserve(n_faces);
    centroids.clear(); centroids.reserve(n_faces);
    
    for (SimpleFace sface : *this) {
        Vec3 v0 = get_vec(sface[0]);
        Vec3 v1 = get_vec(sface[1]);
        Vec3 v2 = get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 n = e1.crossProduct(e2);
        
        areas.push_back(n.norm() * 0.5);
        norms.push_back(n.normalize());
        centroids.push_back((v0 + v1 + v2) / 3.0);
    }

    init_markers(n_faces, TYPES.SURFACE);
    calc_statistics();
}

Point3 TetgenFaces::get_centroid(const int i) const {
    require(i >= 0 && i < static_cast<int>(centroids.size()), "Invalid index: " + d2s(i));
    return centroids[i];
}

Vec3 TetgenFaces::get_norm(const int i) const {
    require(i >= 0 && i < static_cast<int>(norms.size()), "Invalid index: " + d2s(i));
    return norms[i];
}

double TetgenFaces::get_area(const int i) const {
    require(i >= 0 && i < static_cast<int>(areas.size()), "Invalid index: " + d2s(i));
    return areas[i];
}
    
void TetgenFaces::init_statistics() {
    stat.edgemin = DBL_MAX;
    stat.edgemax = 0;
    stat.xmin = stat.ymin = stat.zmin = DBL_MAX;
    stat.xmax = stat.ymax = stat.zmax =-DBL_MAX;
    stat.xbox = stat.ybox = stat.zbox = 0;
}

void TetgenFaces::calc_statistics() {
    init_statistics();

    // loop through all the triangles
    for (SimpleFace face : *this) {
        // read the coordinates of triangle nodes
        Point3 n1 = get_node(face[0]);
        Point3 n2 = get_node(face[1]);
        Point3 n3 = get_node(face[2]);

        // calculate squared lengths of the triangle edges
        const double e1 = n1.distance2(n2);
        const double e2 = n1.distance2(n3);
        const double e3 = n2.distance2(n3);

        // store min and max length
        stat.edgemin = min(stat.edgemin, min(e1, min(e2, e3)));
        stat.edgemax = max(stat.edgemax, max(e1, max(e2, e3)));
    }

    // squared length -> length
    stat.edgemin = sqrt(stat.edgemin);
    stat.edgemax = sqrt(stat.edgemax);

    // calculate the span of centroid coordinates
    for (unsigned int i = 0; i < centroids.size(); ++i) {
        Point3 centroid = get_centroid(i);
        stat.xmax = max(stat.xmax, centroid.x);
        stat.xmin = min(stat.xmin, centroid.x);
        stat.ymax = max(stat.ymax, centroid.y);
        stat.ymin = min(stat.ymin, centroid.y);
        stat.zmax = max(stat.zmax, centroid.z);
        stat.zmin = min(stat.zmin, centroid.z);
    }
    stat.xbox = stat.xmax - stat.xmin;
    stat.ybox = stat.ymax - stat.ymin;
    stat.zbox = stat.zmax - stat.zmin;
}

/* =====================================================================
 *  ========================== TetgenElements =========================
 * ===================================================================== */

void TetgenElements::init_statistics() {
    stat.qmin = stat.edgemin = DBL_MAX;
    stat.qmax = stat.edgemax = -DBL_MAX;
    stat.qmean = 0;
}

void TetgenElements::calc_statistics() {
    init_statistics();

    // loop through all the tetrahedra
    for (SimpleElement elem : *this) {
        // read the coordinates of tetrahedron nodes
        Point3 n1 = get_node(elem[0]);
        Point3 n2 = get_node(elem[1]);
        Point3 n3 = get_node(elem[2]);
        Point3 n4 = get_node(elem[3]);

        // calculate squared lengths of the tetrahedron edges
        const double e1 = n1.distance2(n2);
        const double e2 = n1.distance2(n3);
        const double e3 = n1.distance2(n4);
        const double e4 = n2.distance2(n3);
        const double e5 = n2.distance2(n4);
        const double e6 = n3.distance2(n4);

        // store min and max length
        stat.edgemin = min(stat.edgemin, min(e1, min(e2, min(e3, min(e4, min(e5, e6))))));
        stat.edgemax = max(stat.edgemax, max(e1, max(e2, max(e3, max(e4, max(e5, e6))))));
    }

    const int n_elems = size();

    for (int i = 0; i < n_elems; ++i) {
        double q = get_quality(i);
        stat.qmin = min(stat.qmin, q);
        stat.qmax = max(stat.qmax, q);
        stat.qmean += q;
    }
    stat.qmean /= (1.0*n_elems);

    // squared length -> length
    stat.edgemin = sqrt(stat.edgemin);
    stat.edgemax = sqrt(stat.edgemax);
}

void TetgenElements::transfer(const bool write2read) {
    if (reads == writes) return;
    TetgenCells::transfer(write2read);
    init_statistics();
    if (write2read) {
        if (reads->tetrahedronlist != (int *) NULL)
            delete[] reads->tetrahedronlist;
        reads->tetrahedronlist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            reads->tetrahedronlist[i] = writes->tetrahedronlist[i];
    } else {
        if (writes->tetrahedronlist != (int *) NULL)
            delete[] writes->tetrahedronlist;
        writes->tetrahedronlist = new int[DIM * i_cells];
        for (int i = 0; i < DIM * i_cells; ++i)
            writes->tetrahedronlist[i] = reads->tetrahedronlist[i];
    }
}

void TetgenElements::init(const int N) {
    TetgenCells::init(N);
    init_statistics();
    if (writes->tetrahedronlist != (int *) NULL)
        delete[] writes->tetrahedronlist;
    writes->tetrahedronlist = new int[DIM * N];
}

void TetgenElements::append(const SimpleCell<4> &cell) {
    require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
    int i = DIM * i_cells;
    for (unsigned int node : cell)
        writes->tetrahedronlist[i++] = node;
    i_cells++;
}

double TetgenElements::get_volume(const int i) const {
    SimpleElement elem = get_cell(i);

    // read nodes
    Vec3 n0 = get_vec(elem[0]);
    Vec3 e1 = get_vec(elem[1]);
    Vec3 e2 = get_vec(elem[2]);
    Vec3 e3 = get_vec(elem[3]);

    // calculate edges
    e1 -= n0;
    e2 -= n0;
    e3 -= n0;

    return abs( e1.dotProduct(e2.crossProduct(e3))/6.0 );
}

double TetgenElements::get_quality(const int i) const {
    SimpleElement elem = get_cell(i);

    // read nodes
    Vec3 n0 = get_vec(elem[0]);
    Vec3 e1 = get_vec(elem[1]);
    Vec3 e2 = get_vec(elem[2]);
    Vec3 e3 = get_vec(elem[3]);

    // calculate edges
    Vec3 e4 = e1 - e2;
    Vec3 e5 = e1 - e3;
    Vec3 e6 = e2 - e3;
    e1 -= n0;
    e2 -= n0;
    e3 -= n0;

    // calculate quality
    // no need to divide by 6 top & bot, as they cancel out
    const double top = abs(e1.dotProduct(e2.crossProduct(e3)));
    const double bot = sqrt(e1.norm2()+e2.norm2()+e3.norm2()+e4.norm2()+e5.norm2()+e6.norm2());
    return top / bot;
}

SimpleCell<4> TetgenElements::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    const int I = DIM * i;
    return SimpleElement(reads->tetrahedronlist[I], reads->tetrahedronlist[I+1],
            reads->tetrahedronlist[I+2], reads->tetrahedronlist[I+3]);
}

vector<int> TetgenElements::get_neighbours(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    require(reads->neighborlist, "Query from empty neighbour list!");

    const int I = DIM * i;
    const int* nborlist = reads->neighborlist;
    return vector<int> {nborlist[I+0], nborlist[I+1], nborlist[I+2], nborlist[I+3]};
}

/* =====================================================================
 *  ========================== Quadrangles ===========================
 * ===================================================================== */

void Quadrangles::init(const int N) {
    TetgenCells::init(N);
    quads.clear();
    quads.reserve(N);
}

void Quadrangles::append(const SimpleCell<4> &cell) {
    expect(quads.size() < quads.capacity(), "Allocated size of cells exceeded!");
    quads.push_back(cell);
    i_cells++;
}

SimpleCell<4> Quadrangles::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    return quads[i];
}

/* =====================================================================
 *  ============================ Hexahedra ============================
 * ===================================================================== */

void Hexahedra::init(const int N) {
    TetgenCells::init(N);
    hexs.clear();
    hexs.reserve(N);
}

void Hexahedra::append(const SimpleCell<8> &cell) {
    expect(hexs.size() < hexs.capacity(), "Allocated size of cells exceeded!");
    hexs.push_back(cell);
    i_cells++;
}

SimpleCell<8> Hexahedra::get_cell(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
    return hexs[i];
}

vector<dealii::CellData<3>> Hexahedra::export_vacuum() const {
    const int n_elems = vector_sum(vector_greater(&markers, 0));
    std::vector<dealii::CellData<3>> elems; elems.reserve(n_elems);

    // loop through all the hexahedra and export the ones that are in the vacuum
    for (int i = 0; i < size(); ++i)
        if (get_marker(i) > 0) {
            elems.push_back(dealii::CellData<3>());
            SimpleHex hex = (*this)[i];
            for (int v = 0; v < DIM; ++v)
                elems.back().vertices[v] = hex[v];
        }
    return elems;
}

vector<dealii::CellData<3>> Hexahedra::export_bulk() const {
    const int n_elems = vector_sum(vector_less(&markers, 0));
    std::vector<dealii::CellData<3>> elems; elems.reserve(n_elems);

    // loop through all the hexahedra and export the ones that are in the bulk
    for (int i = 0; i < size(); ++i)
        if (get_marker(i) < 0) {
            elems.push_back(dealii::CellData<3>());
            SimpleHex hex = (*this)[i];
            for (int v = 0; v < DIM; ++v)
                elems.back().vertices[v] = hex[v];
        }
    return elems;
}

} /* namespace femocs */
