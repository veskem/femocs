/*
 * Mesh.cpp
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#include "Mesh.h"
#include <fstream>
#include <float.h>

using namespace std;
namespace femocs {

// Mesh constructor
Mesh::Mesh(const string& mesher) : i_nodes(0), i_edges(0), i_elems(0), i_faces(0) {
    require(mesher == "tetgen", "Unimplemented mesher!");
    this->mesher = mesher;
    stat.Vmin = stat.Vmax = stat.Vmedian = stat.Vaverage = 0;
}

// Mesh destructor
Mesh::~Mesh() {}

// Function to generate simple mesh that consists of one element
const void Mesh::generate_simple() {
    const int n_nodes = n_nodes_per_elem;
    const int n_edges = n_edges_per_elem;
    const int n_faces = n_faces_per_elem;
    const int n_elems = 1;

    init_nodes(n_nodes);
    add_node(Point3(1.0, 0.0, 0.7));
    add_node(Point3(-1.0, 0.0, 0.7));
    add_node(Point3(0.0, 1.0, -0.7));
    add_node(Point3(0.0, -1.0, -0.7));

    init_faces(n_faces);
    add_face(0, 1, 3);
    add_face(1, 2, 3);
    add_face(2, 0, 3);
    add_face(0, 1, 2);

    init_edges(n_edges);
    add_edge(0, 1);
    add_edge(0, 2);
    add_edge(0, 3);
    add_edge(1, 2);
    add_edge(1, 3);
    add_edge(2, 3);

    init_elems(n_elems);
    add_elem(0, 1, 2, 3);

    recalc("rQ");
}

// =================================
// *** GETTERS: ***************

const Vec3 Mesh::get_vec(const int i) const {
    require(get_n_nodes() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
    const int n = n_coordinates * i;
    return Vec3(tetIOout.pointlist[n+0], tetIOout.pointlist[n+1], tetIOout.pointlist[n+2]);
}

const Point3 Mesh::get_node(const int i) const {
    require(get_n_nodes() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
    const int n = n_coordinates * i;
    return Point3(tetIOout.pointlist[n+0], tetIOout.pointlist[n+1], tetIOout.pointlist[n+2]);
}

const SimpleEdge Mesh::get_simpleedge(const int i) const {
    require(get_n_edges() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_edges(), "Invalid index: " + to_string(i));
    const int I = n_nodes_per_edge * i;
    return SimpleEdge(tetIOout.edgelist[I], tetIOout.edgelist[I+1]);
}

const SimpleFace Mesh::get_simpleface(const int i) const {
    require(get_n_faces() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_faces(), "Invalid index: " + to_string(i));
    const int I = n_nodes_per_face * i;
    return SimpleFace(tetIOout.trifacelist[I], tetIOout.trifacelist[I+1], tetIOout.trifacelist[I+2]);
}

const SimpleElement Mesh::get_simpleelem(const int i) const {
    require(get_n_elems() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_elems(), "Invalid index: " + to_string(i));
    const int I = n_nodes_per_elem * i;
    return SimpleElement(tetIOout.tetrahedronlist[I], tetIOout.tetrahedronlist[I+1],
            tetIOout.tetrahedronlist[I+2], tetIOout.tetrahedronlist[I+3]);
}

const Point3 Mesh::get_face_centroid(int i) const {
    require(get_n_faces() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_faces(), "Invalid index: " + to_string(i));

    Point3 verts(0, 0, 0);
    for (int v : get_simpleface(i))
        verts += get_node(v);
    verts /= 1.0*n_nodes_per_face;

    return verts;
}

const Point3 Mesh::get_elem_centroid(int i) const {
    require(get_n_elems() > 0, "Inquiry from empty mesh!");
    require(i >= 0 && i < get_n_elems(), "Invalid index: " + to_string(i));

    Point3 verts;
    for (int v : get_simpleelem(i))
        verts += get_node(v);
    verts *= 1.0/n_nodes_per_elem;

    return verts;
}

const double Mesh::get_area(const int i) const {
    require(i >= 0 && i < get_n_areas(), "Invalid index: " + to_string(i));
    return areas[i];
}

const double Mesh::get_volume(const int i) const {
    require(i >= 0 && i < get_n_volumes(), "Invalid index: " + to_string(i));
    return volumes[i];
}

const double Mesh::get_quality(const int i) const {
    require(i >= 0 && i < get_n_qualities(), "Invalid index: " + to_string(i));
    return qualities[i];
}

const int Mesh::get_nodemarker(const int i) const {
    require(i >= 0 && i < get_n_nodemarkers(), "Invalid index: " + to_string(i));
    return nodemarkers[i];
}

const int Mesh::get_edgemarker(const int i) const {
    require(i >= 0 && i < get_n_edgemarkers(), "Invalid index: " + to_string(i));
    return edgemarkers[i];
}

const int Mesh::get_facemarker(const int i) const {
    require(i >= 0 && i < get_n_facemarkers(), "Invalid index: " + to_string(i));
    return facemarkers[i];
}

const int Mesh::get_elemmarker(const int i) const {
    require(i >= 0 && i < get_n_elemmarkers(), "Invalid index: " + to_string(i));
    return elemmarkers[i];
}

const double* Mesh::get_nodes() {
    return tetIOout.pointlist;
}

const int* Mesh::get_edges() {
    return tetIOout.edgelist;
}

const int* Mesh::get_faces() {
    return tetIOout.trifacelist;
}

const int* Mesh::get_elems() {
    return tetIOout.tetrahedronlist;
}

const vector<int>* Mesh::get_nodemarkers() {
    return &nodemarkers;
}

const vector<int>* Mesh::get_edgemarkers() {
    return &edgemarkers;
}

const vector<int>* Mesh::get_facemarkers() {
    return &facemarkers;
}

const vector<int>* Mesh::get_elemmarkers() {
    return &elemmarkers;
}

const vector<int> Mesh::get_elem_neighbours(const int i) const {
    require(i >= 0 && i < get_n_elems(), "Invalid index: " + to_string(i));

    const int I = n_faces_per_elem * i;
    if (tetIOout.neighborlist)
        return vector<int> {tetIOout.neighborlist[I+0], tetIOout.neighborlist[I+1],
            tetIOout.neighborlist[I+2], tetIOout.neighborlist[I+3]};

    require(false, "Query from empty neighbour list!");
    return vector<int> {i};
}

const int Mesh::get_n_nodes() const {
    return tetIOout.numberofpoints;
}

const int Mesh::get_n_edges() const {
    return tetIOout.numberofedges;
}

const int Mesh::get_n_faces() const {
    return tetIOout.numberoftrifaces;
}

const int Mesh::get_n_elems() const {
    return tetIOout.numberoftetrahedra;
}

const int Mesh::get_n_nodemarkers() const {
    return nodemarkers.size();
}

const int Mesh::get_n_edgemarkers() const {
    return edgemarkers.size();
}

const int Mesh::get_n_facemarkers() const {
    return facemarkers.size();
}

const int Mesh::get_n_elemmarkers() const {
    return elemmarkers.size();
}

const int Mesh::get_n_volumes() const {
    return volumes.size();
}

const int Mesh::get_n_areas() const {
    return areas.size();
}

const int Mesh::get_n_qualities() const {
    return qualities.size();
}

// =================================
// *** SETTERS: ********************
// =================================

const void Mesh::set_nodemarker(const int node, const int m) {
    require(node >= 0 && node < get_n_nodemarkers(), "Invalid index: " + to_string(node));
    nodemarkers[node] = m;
}

const void Mesh::set_edgemarker(const int edge, const int m) {
    require(edge >= 0 && edge < get_n_edgemarkers(), "Invalid index: " + to_string(edge));
    edgemarkers[edge] = m;
}

const void Mesh::set_facemarker(const int face, const int m) {
    require(face >= 0 && face < get_n_facemarkers(), "Invalid index: " + to_string(face));
    facemarkers[face] = m;
}

const void Mesh::set_elemmarker(const int elem, const int m){
    require(elem >= 0 && elem < get_n_elemmarkers(), "Invalid index: " + to_string(elem));
    elemmarkers[elem] = m;
}

// =================================
// *** INITIALIZERS: ***************

const void Mesh::init_nodemarkers(const int N) {
    require(N > 0, "Invalid number of node markers: " + to_string(N));
    nodemarkers.reserve(N);
}

const void Mesh::init_edgemarkers(const int N) {
    require(N > 0, "Invalid number of edge markers: " + to_string(N));
    edgemarkers.reserve(N);
}

const void Mesh::init_facemarkers(const int N) {
    require(N > 0, "Invalid number of face markers: " + to_string(N));
    facemarkers.reserve(N);
}

const void Mesh::init_elemmarkers(const int N) {
    require(N > 0, "Invalid number of element markers: " + to_string(N));
    elemmarkers.reserve(N);
}

const void Mesh::init_nodes(const int N) {
    require(N > 0, "Invalid number of nodes: " + to_string(N));
    i_nodes = 0;
    tetIOin.numberofpoints = N;
    tetIOin.pointlist = new REAL[n_coordinates * N];
}

// NB! Edges are added to tetIOout not to tetIOin.
// That's because inserted edges are intended to be valid immediately without re-calc in tetgen
const void Mesh::init_edges(const int N) {
    require(N > 0, "Invalid number of edges: " + to_string(N));
    i_edges = 0;
    tetIOout.numberofedges = N;
    tetIOout.edgelist = new int[n_nodes_per_edge * N];
}

const void Mesh::init_faces(const int N) {
    require(N > 0, "Invalid number of faces: " + to_string(N));
    i_faces = 0;
    tetIOout.numberoftrifaces = N;
    tetIOout.trifacelist = new int[n_nodes_per_face * N];
}

const void Mesh::init_elems(const int N) {
    require(N > 0, "Invalid number of elements: " + to_string(N));
    i_elems = 0;
    tetIOin.numberoftetrahedra = N;
    tetIOin.tetrahedronlist = new int[n_nodes_per_elem * N];
}

// =================================
// *** ADDERS: ***************

const void Mesh::add_node(const double x, const double y, const double z) {
    add_node(Point3(x, y, z));
}

const void Mesh::add_edge(const int n1, const int n2) {
    add_edge(SimpleEdge(n1, n2));
}

const void Mesh::add_face(const int n1, const int n2, const int n3) {
    add_face(SimpleFace(n1, n2, n3));
}

const void Mesh::add_elem(const int n1, const int n2, const int n3, const int n4) {
    add_elem(SimpleElement(n1, n2, n3, n4));
}

const void Mesh::add_node(const Point3 &point) {
    require(i_nodes < tetIOin.numberofpoints, "Allocated size of nodes exceeded!");
    int i = n_coordinates * i_nodes;
    for (double node : point)
        tetIOin.pointlist[i++] = node;
    i_nodes++;
}

// NB! Edges are added to tetIOout not to tetIOin.
// That's because inserted edges are intended to be valid immediately without re-calc in tetgen
const void Mesh::add_edge(const SimpleEdge& edge) {
    require(i_edges < tetIOout.numberofedges, "Allocated size of edges exceeded!");
    require(edge.n1 >= 0 && edge.n2 >= 0, "Invalid edge: " + edge.to_str());
    int i = n_nodes_per_edge * i_edges;
    for (int node : edge)
        tetIOout.edgelist[i++] = node;
    i_edges++;
}

const void Mesh::add_face(const SimpleFace& face) {
    require(i_faces < tetIOout.numberoftrifaces, "Allocated size of faces exceeded!");
    require(face.n1 >= 0 && face.n2 >= 0 && face.n3 >= 0, "Invalid face: " + face.to_str());
    int i = n_nodes_per_face * i_faces;
    for (int node : face)
        tetIOout.trifacelist[i++] = node;
    i_faces++;
}

const void Mesh::add_elem(const SimpleElement& el) {
    require(i_elems < tetIOin.numberoftetrahedra, "Allocated size of elements exceeded!");
    require(el.n1 >= 0 && el.n2 >= 0 && el.n3 >= 0 && el.n4 >= 0, "Invalid element: " + el.to_str());
    int i = n_nodes_per_elem * i_elems;
    for (int node : el)
        tetIOin.tetrahedronlist[i++] = node;
    i_elems++;
}

const void Mesh::add_nodemarker(const int m) {
    expect(get_n_nodemarkers() < nodemarkers.capacity(), "Allocated size of nodemarkers exceeded!");
    nodemarkers.push_back(m);
}

const void Mesh::add_edgemarker(const int m) {
    expect(get_n_edgemarkers() < edgemarkers.capacity(), "Allocated size of edgemarkers exceeded!");
    edgemarkers.push_back(m);
}

const void Mesh::add_facemarker(const int m) {
    expect(get_n_facemarkers() < facemarkers.capacity(), "Allocated size of facemarkers exceeded!");
    facemarkers.push_back(m);
}

const void Mesh::add_elemmarker(const int m) {
    expect(get_n_elemmarkers() < elemmarkers.capacity(), "Allocated size of elemmarkers exceeded!");
    elemmarkers.push_back(m);
}

// =================================
// *** REPLICATORS: ***************

const void Mesh::copy_statistics(Mesh* mesh) {
    stat.Vmin = mesh->stat.Vmin;
    stat.Vmax = mesh->stat.Vmax;
    stat.Vaverage = mesh->stat.Vaverage;
    stat.Vmedian = mesh->stat.Vmedian;
    stat.n_bulk = mesh->stat.n_bulk;
    stat.n_surface = mesh->stat.n_surface;
    stat.n_vacuum = mesh->stat.n_vacuum;
}

const void Mesh::copy_nodes(Mesh* mesh) {
    int n_nodes = mesh->get_n_nodes();
    const double* nodes = mesh->get_nodes();

    init_nodes(n_nodes);
    for (int i = 0; i < n_coordinates * n_nodes; ++i)
        tetIOin.pointlist[i] = nodes[i];
    i_nodes = n_nodes;
}

const void Mesh::copy_nodes(Mesh* mesh, const vector<bool> &mask) {
    init_nodes(vector_sum(mask));
    for (int i = 0; i < mesh->get_n_nodes(); ++i)
        if (mask[i])
            add_node(mesh->get_node(i));
}

const void Mesh::copy_edges(Mesh* mesh) {
    const int n_edges = mesh->get_n_edges();
    init_edges(n_edges);
    for (int i = 0; i < n_edges; ++i)
        add_edge(mesh->get_simpleedge(i));
}

const void Mesh::copy_edges(const vector<SimpleEdge> &edges, const vector<bool> &mask) {
    const int n_edges = edges.size();
    require(n_edges == mask.size(), "Incompatible vector sizes!");

    init_edges(vector_sum(mask));
    for (int i = 0; i < n_edges; ++i)
        if (mask[i])
            add_edge(edges[i]);
}

const void Mesh::copy_faces(Mesh* mesh) {
    int n_faces = mesh->get_n_faces();
    const int* faces = mesh->get_faces();

    init_nodes(n_faces);
    for (int i = 0; i < n_nodes_per_face * n_faces; ++i)
        tetIOin.trifacelist[i] = faces[i];
    i_faces = n_faces;
}

const void Mesh::copy_elems(Mesh* mesh) {
    int n_elems = mesh->get_n_elems();
    const int* elems = mesh->get_elems();

    init_nodes(n_elems);
    for (int i = 0; i < n_nodes_per_elem * n_elems; ++i)
        tetIOin.tetrahedronlist[i] = elems[i];
    i_elems = n_elems;
}

const void Mesh::copy_nodemarkers(Mesh* mesh) {
    int N = mesh->get_n_nodemarkers();
    for (int i = 0; i < N; ++i)
        nodemarkers.push_back(mesh->get_nodemarker(i));
}

const void Mesh::copy_edgemarkers(Mesh* mesh) {
    int N = mesh->get_n_edgemarkers();
    for (int i = 0; i < N; ++i)
        edgemarkers.push_back(mesh->get_edgemarker(i));
}

const void Mesh::copy_facemarkers(Mesh* mesh) {
    int N = mesh->get_n_facemarkers();
    for (int i = 0; i < N; ++i)
        facemarkers.push_back(mesh->get_facemarker(i));
}

const void Mesh::copy_elemmarkers(Mesh* mesh) {
    int N = mesh->get_n_elemmarkers();
    for (int i = 0; i < N; ++i)
        elemmarkers.push_back(mesh->get_elemmarker(i));
}

// =================================
// *** VARIA: ***************

const Medium Mesh::to_medium() const {
    int n_atoms = get_n_nodes();
    Medium medium;
    medium.reserve(n_atoms);

    for(int i = 0; i < n_atoms; ++i)
        medium.add_atom( Atom(i, get_node(i), 0) );

    medium.calc_statistics();
    return medium;
}

const double Mesh::calc_volume(const int i) {
    require(i >= 0 && i < get_n_elems(), "Invalid index: " + to_string(i));

    SimpleElement elem = get_simpleelem(i);
    Vec3 base = get_vec(elem[0]);
    Vec3 edge1 = base - get_vec(elem[1]);
    Vec3 edge2 = base - get_vec(elem[2]);
    Vec3 edge3 = base - get_vec(elem[3]);

    double volume = edge1.dotProduct(edge2.crossProduct(edge3));
    return fabs(volume) / 6.0;
}

const double Mesh::calc_area(const int i) {
    require(i >= 0 && i < get_n_faces(), "Invalid index: " + to_string(i));

    SimpleFace face = get_simpleface(i);
    Vec3 base = get_vec(face[0]);
    Vec3 edge1 = base - get_vec(face[1]);
    Vec3 edge2 = base - get_vec(face[2]);

    double area = edge1.crossProduct(edge2).length();
    return fabs(area) / 2.0;
}

const void Mesh::calc_volumes() {
    int elem;
    int n_elems = get_n_elems();

    volumes.reserve(n_elems);

    // Loop through the elements
    for (elem = 0; elem < n_elems; ++elem)
        volumes.push_back(calc_volume(elem));
}

const void Mesh::calc_areas() {
    int face;
    int n_faces = get_n_faces();

    areas.reserve(n_faces);

    // Loop through the faces
    for (face = 0; face < n_faces; ++face)
        areas.push_back(calc_area(face));
}

// Return determinant of matrix composed of n1, n2 and n3
const double Mesh::determinant(const Vec3 &n1, const Vec3 &n2, const Vec3 &n3) {
    return n1.x * (n2.y*n3.z - n2.z*n3.y)
            - n1.y * (n1.x*n3.z - n1.z*n3.x)
            + n1.z * (n1.x*n2.y - n1.y*n2.x);
}

// Return determinant of matrix composed of n1, n2, n3 and n4
double Mesh::determinant(const double* n1, const double* n2, const double* n3, const double* n4, bool ones) {
    double det1 = determinant(Vec3(n1[1],n1[2],n1[3]), Vec3(n2[1],n2[2],n2[3]), Vec3(n3[1],n3[2],n3[3]));
    double det2 = determinant(Vec3(n1[0],n1[2],n1[3]), Vec3(n2[0],n2[2],n2[3]), Vec3(n3[0],n3[2],n3[3]));
    double det3 = determinant(Vec3(n1[0],n1[1],n1[3]), Vec3(n2[0],n2[1],n2[3]), Vec3(n3[0],n3[1],n3[3]));
    double det4 = determinant(Vec3(n1[0],n1[1],n1[2]), Vec3(n2[0],n2[1],n2[2]), Vec3(n3[0],n3[1],n3[2]));

    if(ones) return det1 - det2 + det3 - det4;
    else return n4[0] * det1 - n4[1] * det2 + n4[2] * det3 - n4[3] * det4;
}

// circumsphere of tetrahedron: http://mathworld.wolfram.com/Circumsphere.html
const void Mesh::calc_qualities_byelem() {
    const int n_elems = get_n_elems();
    qualities.resize(n_elems);

    for(int el = 0; el < n_elems; ++el) {
        SimpleElement elem = get_simpleelem(el);
        Vec3 node1 = get_vec(elem.n1);
        Vec3 node2 = get_vec(elem.n2);
        Vec3 node3 = get_vec(elem.n3);
        Vec3 node4 = get_vec(elem.n4);

        // Tetrahedron edge lengths
        double e1 = (node1 - node2).length();
        double e2 = (node1 - node3).length();
        double e3 = (node1 - node4).length();
        double e4 = (node2 - node3).length();
        double e5 = (node2 - node4).length();
        double e6 = (node3 - node4).length();

        // The length of the tetrahedron shortest edge
        double min_edge = min(e1, min(e2, min(e3, min(e4, min(e5, e6)))));

        double norm[] = {node1.length(), node2.length(), node3.length(), node4.length()};
        double xx[] = {node1.x, node2.x, node3.x, node4.x};
        double yy[] = {node1.y, node2.y, node3.y, node4.y};
        double zz[] = {node1.z, node2.z, node3.z, node4.z};
        double one[] = {1.0, 1.0, 1.0, 1.0};

        double a  = determinant(xx, yy, zz, one, true);
        double Dx = determinant(norm, yy, zz, one, true);
        double Dy = -1.0*determinant(norm, xx, zz, one, true);
        double Dz = determinant(norm, xx, yy, one, true);
        double c  = determinant(norm, xx, yy, zz, false);

        // The radius of tetrahedron circumsphere
        double r = sqrt(Dx*Dx + Dy*Dy + Dz*Dz - 4*a*c) / (2*fabs(a));
        qualities[el] = (r / min_edge);
    }
}

// Quality of elem = max(quality of face) = max(circumradius / min(length of edge))
// circumradius of triangle: http://mathworld.wolfram.com/Circumradius.html
const void Mesh::calc_qualities_byface() {
    const int n_elems = get_n_elems();
    qualities.resize(n_elems);

    for(int el = 0; el < n_elems; ++el) {
        SimpleElement elem = get_simpleelem(el);
        Vec3 node1 = get_vec(elem.n1);
        Vec3 node2 = get_vec(elem.n2);
        Vec3 node3 = get_vec(elem.n3);
        Vec3 node4 = get_vec(elem.n4);

        double q1 = calc_face_quality(node1, node2, node3);
        double q2 = calc_face_quality(node1, node2, node4);
        double q3 = calc_face_quality(node1, node3, node4);
        double q4 = calc_face_quality(node2, node3, node4);
        qualities[el] = (max(q1, max(q2, max(q3, q4))));
    }
}

const double Mesh::calc_face_quality(const Vec3 &node1, const Vec3 &node2, const Vec3 &node3) {
    double a = (node1 - node2).length();
    double b = (node1 - node3).length();
    double c = (node2 - node3).length();
    double D = (a + b + c) * (b + c - a) * (c + a - b) * (a + b - c);

    double min_edge = min(a, min(b, c));
    double R = a * b * c / sqrt(D);

    return R / min_edge;
}

const void Mesh::init_statistics() {
    stat.n_bulk = stat.n_surface = stat.n_vacuum = 0;
    stat.Vmin = stat.Vmax = stat.Vaverage = stat.Vmedian = 0;
    stat.xmin = stat.ymin = stat.zmin = DBL_MAX;
    stat.xmax = stat.ymax = stat.zmax = DBL_MIN;
}

const void Mesh::calc_statistics() {
    init_statistics();
    size_t n_nodes = get_n_nodes();

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
}

const void Mesh::calc_statistics(const int i) {
    size_t n_volumes = get_n_volumes();
    size_t n_markers = get_n_nodemarkers();

    calc_statistics();

    // Find the number of nodes in various regions
    if (n_markers > 0)
        // Loop through all the node markers
        for (int marker : nodemarkers) {
            if (marker == TYPES.BULK)
                stat.n_bulk++;
            else if (marker == TYPES.VACUUM)
                stat.n_vacuum++;
            else if (marker == TYPES.SURFACE)
                stat.n_surface++;
        }

    // Calculate the statistics about the volumes of elements
    if (n_volumes > 0) {
        // Make a copy of the volumes vector
        vector<double> tempvec;
        tempvec.reserve(n_volumes);
        copy(volumes.begin(), volumes.end(), back_inserter(tempvec));

        sort(tempvec.begin(), tempvec.end());

        // Store min, max and average volume
        stat.Vmin = tempvec[0];
        stat.Vmax = tempvec[n_volumes - 1];
        stat.Vaverage = vector_sum(tempvec) / n_volumes;

        // Store the median volume
        if (n_volumes % 2 == 0)
            stat.Vmedian = (tempvec[n_volumes / 2 - 1] + tempvec[n_volumes / 2]) / 2;
        else
            stat.Vmedian = tempvec[n_volumes / 2];
    }
}

const void Mesh::recalc() {
    i_nodes = tetIOin.numberofpoints;
    i_faces = tetIOin.numberoftrifaces;
    i_elems = tetIOin.numberoftetrahedra;

    tetIOout.numberofpoints = i_nodes;
    tetIOout.numberoftrifaces = i_faces;
    tetIOout.numberoftetrahedra = i_elems;

    tetIOout.pointlist = new REAL[n_coordinates * i_nodes];
    tetIOout.trifacelist = new int[n_nodes_per_face * i_faces];
    tetIOout.tetrahedronlist = new int[n_nodes_per_elem * i_elems];

    for (int i = 0; i < n_coordinates * i_nodes; ++i)
        tetIOout.pointlist[i] = tetIOin.pointlist[i];

    for (int i = 0; i < n_nodes_per_face * i_faces; ++i)
        tetIOout.trifacelist[i] = tetIOin.trifacelist[i];

    for (int i = 0; i < n_nodes_per_elem * i_elems; ++i)
        tetIOout.tetrahedronlist[i] = tetIOin.tetrahedronlist[i];
}

// Function to perform tetgen calculation on input and write_tetgen data
const void Mesh::recalc(const string& cmd) {
    tetrahedralize(const_cast<char*>(cmd.c_str()), &tetIOin, &tetIOout);

    i_nodes = tetIOout.numberofpoints;
    i_faces = tetIOout.numberoftrifaces;
    i_elems = tetIOout.numberoftetrahedra;
}

// Function to perform tetgen calculation on input and write_tetgen data
const void Mesh::recalc(const string& cmd1, const string& cmd2) {
    tetgenio tetIOtemp;

    tetrahedralize(const_cast<char*>(cmd1.c_str()), &tetIOin, &tetIOtemp);
    tetrahedralize(const_cast<char*>(cmd2.c_str()), &tetIOtemp, &tetIOout);

    i_nodes = tetIOout.numberofpoints;
    i_faces = tetIOout.numberoftrifaces;
    i_elems = tetIOout.numberoftetrahedra;
}

// Write tetgen mesh into files with its internal functions
const void Mesh::write_tetgen(const string file_name) {
    const string cmd = "Q";
    tetgenbehavior tetgenbeh;

    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    for (int i = 0; i < file_name.size(); ++i)
        tetgenbeh.outfilename[i] = file_name[i];

    tetrahedralize(&tetgenbeh, &tetIOout, NULL);
}

// Function to write_tetgen mesh in .vtk format
const void Mesh::write_vtk(const string file_name, const int n_nodes, const int n_cells,
        const int n_markers, const REAL* nodes, const int* cells, const vector<int>* markers,
        const int celltype, const int n_nodes_in_cell) {

    int i, j;
    char file_name_char[1024];
    strcpy(file_name_char, file_name.c_str());

    expect(n_nodes > 0, "Zero nodes detected!");

    FILE *out_file;
    out_file = fopen(file_name_char, "w");
    require(out_file != (FILE*) NULL, "Can't open a file " + file_name);

    fprintf(out_file, "# vtk DataFile Version 3.0\n");
    fprintf(out_file, "# Unstructured grid\n");
    fprintf(out_file, "ASCII\n"); // another option is BINARY
    fprintf(out_file, "DATASET UNSTRUCTURED_GRID\n\n");

    // Output the nodes
    if (n_nodes > 0) {
        fprintf(out_file, "POINTS %d double\n", n_nodes);
        for (i = 0; i < n_coordinates * n_nodes; i += n_coordinates)
            fprintf(out_file, "%.8g %.8g %.8g\n", nodes[i + 0], nodes[i + 1], nodes[i + 2]);
        fprintf(out_file, "\n");
    }

    // Output the cells (tetrahedra, triangles or vertices)
    if (n_cells > 0) {
        fprintf(out_file, "CELLS %d %d\n", n_cells, n_cells * (n_nodes_in_cell + 1));
        for (i = 0; i < n_nodes_in_cell * n_cells; i += n_nodes_in_cell) {
            fprintf(out_file, "%d ", n_nodes_in_cell);
            for (j = 0; j < n_nodes_in_cell; ++j)
                fprintf(out_file, "%d ", cells[i + j]);
            fprintf(out_file, "\n");
        }
        fprintf(out_file, "\n");
    }

    // Output the types of cells, 10=tetrahedron, 5=triangle, 1=vertex
    if (n_cells > 0) {
        fprintf(out_file, "CELL_TYPES %d\n", n_cells);
        for (i = 0; i < n_cells; ++i)
            fprintf(out_file, "%d \n", celltype);
        fprintf(out_file, "\n\n");
    }

    // Output cell attributes
    if ((n_markers > 0) && (n_markers == n_cells)) {
        fprintf(out_file, "CELL_DATA %d\n", n_markers);
        fprintf(out_file, "SCALARS Cell_markers int\n");
        fprintf(out_file, "LOOKUP_TABLE default\n");
        for (i = 0; i < n_markers; ++i)
            fprintf(out_file, "%d\n", (*markers)[i]);
        fprintf(out_file, "\n");
    }

    fclose(out_file);
}

// Function to write_tetgen nodes in .xyz format
const void Mesh::write_xyz(const string file_name) {
    const int n_nodes = get_n_nodes();
    const int n_markers = get_n_nodemarkers();

    expect(n_nodes > 0, "Zero nodes detected!");

    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << get_n_nodes() << "\n";
    out_file << "Mesh nodes: id x y z marker\n";

    for (int i = 0; i < n_nodes; ++i) {
        out_file << i << " ";
        out_file << get_node(i) << " ";

        if(n_nodes == n_markers)
            out_file << get_nodemarker(i) << " ";
        else
            out_file << 0 << " ";

        out_file << endl;
    }

    out_file.close();
}

// Function to write edges to .vtk file
const void Mesh::write_edges(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unknown file type: " + file_type);

    const int celltype = 3; // 1-vertex, 3-line, 5-triangle, 10-tetrahedron

    const int n_nodes = get_n_nodes();
    const int n_edges = get_n_edges();
    const int n_markers = get_n_edgemarkers();
    const REAL* nodes = get_nodes();          // pointer to nodes
    const int* edges = get_edges();           // pointer to edge nodes
    vector<int>* cellmarkers = &edgemarkers;  // pointer to edge markers

    write_vtk(file_name, n_nodes, n_edges, n_markers, nodes, edges, cellmarkers, celltype,
            n_nodes_per_edge);
}

// Function to write faces to .vtk file
const void Mesh::write_faces(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unknown file type: " + file_type);

    const int celltype = 5; // 1-vertex, 5-triangle, 10-tetrahedron

    const int n_nodes = get_n_nodes();
    const int n_faces = get_n_faces();
    const int n_markers = get_n_facemarkers();
    const REAL* nodes = get_nodes();          // pointer to nodes
    const int* faces = get_faces();           // pointer to face nodes
    vector<int>* cellmarkers = &facemarkers;    // pointer to face markers

    write_vtk(file_name, n_nodes, n_faces, n_markers, nodes, faces, cellmarkers, celltype,
            n_nodes_per_face);
}

// Function to write elements to .vtk file
const void Mesh::write_elems(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unknown file type: " + file_type);

    const int celltype = 10; // 1-vertex, 5-triangle, 10-tetrahedron

    const int n_nodes = get_n_nodes();
    const int n_elems = get_n_elems();
    const int n_markers = get_n_elemmarkers();
    const REAL* nodes = get_nodes();          // pointer to nodes
    const int* elems = get_elems();           // pointer to faces nodes
    vector<int>* cellmarkers = &elemmarkers;    // pointer to element markers

    write_vtk(file_name, n_nodes, n_elems, n_markers, nodes, elems, cellmarkers, celltype,
            n_nodes_per_elem);
}

// Function to write nodes to .xyz or .vtk file
const void Mesh::write_nodes(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "xyz" || file_type == "vtk", "Unknown file type: " + file_type);

    if (file_type == "xyz")
        write_xyz(file_name);

    else {
        const int celltype = 1;     // 1-vertex, 3-line, 5-triangle, 10-tetrahedron
        const int n_nodes = get_n_nodes();
        const int n_markers = get_n_nodemarkers();
        const REAL* nodes = get_nodes();           // pointer to nodes
        vector<int>* markers = &nodemarkers; // pointer to node markers

        int node_indx[n_nodes];
        for(int i = 0; i < n_nodes; ++i)
            node_indx[i] = i;

        write_vtk(file_name, n_nodes, n_nodes, n_markers, nodes, node_indx, markers, celltype, 1);
    }

}

// =================================
} /* namespace femocs */
