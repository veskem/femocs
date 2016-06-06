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
Mesh::Mesh(const string mesher) {
    require(mesher == "tetgen", "Unimplemented mesher!");

    this->mesher = mesher;
    i_nodes = 0;
    i_elems = 0;
    i_faces = 0;

    stat.Vmin = stat.Vmax = stat.Vmedian = stat.Vaverage = 0;
}

// Mesh destructor
Mesh::~Mesh() {
}

// Function to generate simple mesh that consists of one element
const void Mesh::generate_simple() {
    const int n_nodes = 4;
    const int n_faces = 4;
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

    init_elems(n_elems);
    add_elem(0, 1, 2, 3);

    recalc("rQ");
}

// =================================
// *** GETTERS: ***************

const Vec3 Mesh::get_vec(const int i) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    const int n = n_coordinates * i;
    return Vec3(tetIOout.pointlist[n+0], tetIOout.pointlist[n+1], tetIOout.pointlist[n+2]);
}

const Point3 Mesh::get_node(const int i) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    const int n = n_coordinates * i;
    return Point3(tetIOout.pointlist[n+0], tetIOout.pointlist[n+1], tetIOout.pointlist[n+2]);
}

const SimpleFace Mesh::get_simpleface(const int i) {
    require(i >= 0 && i < get_n_faces(), "Invalid index!");
    const int I = n_nodes_per_face * i;
    return SimpleFace(tetIOout.trifacelist[I], tetIOout.trifacelist[I+1], tetIOout.trifacelist[I+2]);
}

const SimpleElement Mesh::get_simpleelem(const int i) {
    require(i >= 0 && i < get_n_elems(), "Invalid index!");
    const int I = n_nodes_per_elem * i;
    return SimpleElement(tetIOout.tetrahedronlist[I], tetIOout.tetrahedronlist[I+1], tetIOout.tetrahedronlist[I+2], tetIOout.tetrahedronlist[I+3]);
}

const Point3 Mesh::get_face_centre(int i) {
    require(i >= 0 && i < get_n_faces(), "Invalid index!");

    Point3 verts(0, 0, 0);
    SimpleFace face = get_simpleface(i);
    for (int v = 0; v < n_nodes_per_face; ++v)
        verts += get_node(face[v]);
    verts /= 1.0*n_nodes_per_face;

    return verts;
}

const Point3 Mesh::get_elem_centre(int i) {
    require(i >= 0 && i < get_n_elems(), "Invalid index!");

    Point3 verts;
    SimpleElement elem = get_simpleelem(i);
    for (int v = 0; v < n_nodes_per_elem; ++v)
        verts += get_node(elem[v]);
    verts /= 1.0*n_nodes_per_elem;

    return verts;
}

const double Mesh::get_area(const int i) {
    require(i >= 0 && i < get_n_areas(), "Invalid index!");
    return areas[i];
}

const double Mesh::get_volume(const int i) {
    require(i >= 0 && i < get_n_volumes(), "Invalid index!");
    return volumes[i];
}

const int Mesh::get_nodemarker(const int i) {
    require(i >= 0 && i < get_n_nodemarkers(), "Invalid index!");
    return nodemarkers[i];
}

const int Mesh::get_facemarker(const int i) {
    require(i >= 0 && i < get_n_facemarkers(), "Invalid index!");
    return facemarkers[i];
}

const int Mesh::get_elemmarker(const int i) {
    require(i >= 0 && i < get_n_elemmarkers(), "Invalid index!");
    return elemmarkers[i];
}

const double* Mesh::get_nodes() {
    return tetIOout.pointlist;
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

const vector<int>* Mesh::get_facemarkers() {
    return &facemarkers;
}

const vector<int>* Mesh::get_elemmarkers() {
    return &elemmarkers;
}

const int Mesh::get_n_nodes() {
    return i_nodes;
    //return tetIO.numberofpoints;
}

const int Mesh::get_n_faces() {
    return i_faces;
    //return tetIO.numberoftrifaces;
}

const int Mesh::get_n_elems() {
    return i_elems;
    //return tetIO.numberoftetrahedra;
}

const int Mesh::get_n_nodemarkers() {
    return nodemarkers.size();
}

const int Mesh::get_n_facemarkers() {
    return facemarkers.size();
}

const int Mesh::get_n_elemmarkers() {
    return elemmarkers.size();
}

const int Mesh::get_n_volumes() {
    return volumes.size();
}

const int Mesh::get_n_areas() {
    return areas.size();
}
// =================================
// *** SETTERS: ********************
// =================================

const void Mesh::set_nodemarker(const int node, const int m) {
    require(node >= 0 && node < get_n_nodemarkers(), "Invalid index!");
    nodemarkers[node] = m;
}
const void Mesh::set_facemarker(const int face, const int m) {
    require(face >= 0 && face < get_n_facemarkers(), "Invalid index!");
    facemarkers[face] = m;
}
const void Mesh::set_elemmarker(const int elem, const int m){
    require(elem >= 0 && elem < get_n_elemmarkers(), "Invalid index!");
    elemmarkers[elem] = m;
}

// =================================
// *** INITIALIZERS: ***************

const void Mesh::init_nodemarkers(const int N) {
    require(N > 0, "Invalid number of node markers!");
    nodemarkers.reserve(N);
}

const void Mesh::init_facemarkers(const int N) {
    require(N > 0, "Invalid number of face markers!");
    facemarkers.reserve(N);
}

const void Mesh::init_elemmarkers(const int N) {
    require(N > 0, "Invalid number of element markers!");
    elemmarkers.reserve(N);
}

const void Mesh::init_nodes(const int N) {
    require(N > 0, "Invalid number of nodes!");
    i_nodes = 0;
    tetIOin.numberofpoints = N;
    tetIOin.pointlist = new REAL[3 * N];
}

const void Mesh::init_faces(const int N) {
    require(N > 0, "Invalid number of faces!");
    i_faces = 0;
    tetIOin.numberoftrifaces = N;
    tetIOin.trifacelist = new int[3 * N];
}

const void Mesh::init_elems(const int N) {
    require(N > 0, "Invalid number of elements!");
    i_elems = 0;
    tetIOin.numberoftetrahedra = N;
    tetIOin.tetrahedronlist = new int[4 * N];
}

// =================================
// *** ADDERS: ***************

const void Mesh::add_node(const double n1, const double n2, const double n3) {
    require(get_n_nodes() < tetIOin.numberofpoints, "Allocated size of nodes exceeded!");
    int i = 3 * i_nodes;
    tetIOout.pointlist[i] = n1;
    tetIOout.pointlist[++i] = n2;
    tetIOout.pointlist[++i] = n3;
    i_nodes++;
}

const void Mesh::add_face(const int f1, const int f2, const int f3) {
    require(get_n_faces() < tetIOin.numberoftrifaces, "Allocated size of faces exceeded!");
    int i = 3 * i_faces;
    tetIOin.trifacelist[i] = f1;
    tetIOin.trifacelist[++i] = f2;
    tetIOin.trifacelist[++i] = f3;
    i_faces++;
}

const void Mesh::add_elem(const int e1, const int e2, const int e3, const int e4) {
    require(get_n_elems() < tetIOin.numberoftetrahedra, "Allocated size of elements exceeded!");
    int i = 4 * i_elems;
    tetIOin.tetrahedronlist[i] = e1;
    tetIOin.tetrahedronlist[++i] = e2;
    tetIOin.tetrahedronlist[++i] = e3;
    tetIOin.tetrahedronlist[++i] = e4;
    i_elems++;
}

const void Mesh::add_node(const Point3 &point) {
    require(get_n_nodes() < tetIOin.numberofpoints, "Allocated size of nodes exceeded!");
    int cntr = 0;
    for (int i = n_coordinates * i_nodes; i < n_coordinates * (i_nodes + 1); ++i)
        tetIOin.pointlist[i] = point[cntr++];
    i_nodes++;
}

const void Mesh::add_face(const SimpleFace& face) {
    require(get_n_faces() < tetIOin.numberoftrifaces, "Allocated size of faces exceeded!");
    int cntr = 0;
    for (int i = n_nodes_per_face * i_faces; i < n_nodes_per_face * (i_faces + 1); ++i)
        tetIOin.tetrahedronlist[i] = face[cntr++];
    i_faces++;
}

const void Mesh::add_elem(const SimpleElement& elem) {
    require(get_n_elems() < tetIOin.numberoftetrahedra, "Allocated size of elements exceeded!");
    int cntr = 0;
    for (int i = n_nodes_per_elem * i_elems; i < n_nodes_per_elem * (i_elems+1); ++i)
        tetIOin.tetrahedronlist[i] = elem[cntr++];
    i_elems++;
}

const void Mesh::add_nodemarker(const int m) {
    expect(get_n_nodemarkers() < nodemarkers.capacity(), "Allocated size of nodemarkers exceeded!");
    nodemarkers.push_back(m);
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

const Medium Mesh::to_medium() {
    int n_atoms = get_n_nodes();
    Medium medium;
    medium.reserve(n_atoms);

    for(int i = 0; i < n_atoms; ++i)
        medium.add_atom(i, get_node(i), 0);

    medium.calc_statistics();
    return medium;
}

const double Mesh::calc_volume(const int i) {
    require(i >= 0 && i < get_n_elems(), "Invalid index!");

    SimpleElement elem = get_simpleelem(i);
    Vec3 base = get_vec(elem[0]);
    Vec3 edge1 = base - get_vec(elem[1]);
    Vec3 edge2 = base - get_vec(elem[2]);
    Vec3 edge3 = base - get_vec(elem[3]);

    double volume = edge1.dotProduct(edge2.crossProduct(edge3));
    return fabs(volume) / 6.0;
}

const double Mesh::calc_area(const int i) {
    require(i >= 0 && i < get_n_faces(), "Invalid index!");

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

const void Mesh::init_statistics() {
    stat.n_bulk = stat.n_surface = stat.n_vacuum = 0;
    stat.Vmin = stat.Vmax = stat.Vaverage = stat.Vmedian = 0;
    stat.xmin = stat.ymin = stat.zmin = DBL_MAX;
    stat.xmax = stat.ymax = stat.zmax = DBL_MIN;
}

const void Mesh::calc_statistics(const AtomReader::Types *types) {
    init_statistics();
    size_t n_volumes = get_n_volumes();
    size_t n_markers = get_n_nodemarkers();
    size_t n_nodes = get_n_nodes();

    // Find the min and max coordinates of all nodes
    if (n_nodes > 0)
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
        for (int marker : nodemarkers) {
            if (marker == types->type_bulk)
                stat.n_bulk++;
            else if (marker == types->type_vacuum)
                stat.n_vacuum++;
            else if (marker == types->type_surf)
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

// Function to perform tetgen calculation on input and write_tetgen data
const void Mesh::recalc(const string cmd) {
    tetrahedralize(const_cast<char*>(cmd.c_str()), &tetIOin, &tetIOout);

    i_nodes = tetIOout.numberofpoints;
    i_faces = tetIOout.numberoftrifaces;
    i_elems = tetIOout.numberoftetrahedra;
}

// Function to perform tetgen calculation on input and write_tetgen data
const void Mesh::double_recalc(const string cmd1, const string cmd2) {
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
        for (i = 0; i < 3 * n_nodes; i += 3)
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
    if (n_markers > 0) {
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
    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    const int n_nodes = get_n_nodes();
    const int n_markers = get_n_nodemarkers();

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

// Function to write faces to .vtk file
const void Mesh::write_faces(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string file_type = get_file_type(file_name);
    require(file_type == "vtk", "Unknown file type!");

    const int celltype = 5; // 1-vertex, 5-triangle, 10-tetrahedron

    int n_nodes = get_n_nodes();
    int n_faces = get_n_faces();
    int n_markers = get_n_facemarkers();
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
    require(file_type == "vtk", "Unknown file type!");

    const int celltype = 10; // 1-vertex, 5-triangle, 10-tetrahedron

    int n_nodes = get_n_nodes();
    int n_elems = get_n_elems();
    int n_markers = get_n_elemmarkers();
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
    require(file_type == "xyz" || file_type == "vtk", "Unknown file type!");

    if (file_type == "xyz")
        write_xyz(file_name);

    else {
        const int celltype = 1;     // 1-vertex, 5-triangle, 10-tetrahedron

        int n_nodes = get_n_nodes();
        int n_markers = get_n_nodemarkers();
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
