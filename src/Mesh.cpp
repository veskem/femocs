/*
 * Mesh.cpp
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#include "Mesh.h"

#include <fstream>
#include <algorithm>

using namespace std;
namespace femocs {

// Mesh constructor
Mesh::Mesh(const string mesher) {
    require(mesher == "tetgen", "Unimplemented mesher!");

    i_nodes = 0;
    i_elems = 0;
    i_faces = 0;

    stat.Vmin = stat.Vmax = stat.Vmedian = stat.Vaverage = 0;
}

// Mesh destructor
Mesh::~Mesh() {
}

// =================================
// *** GETTERS: ***************

const Point3d Mesh::get_point(const int i) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    const int n = n_coordinates * i;
    return Point3d(tetIO.pointlist[n+0], tetIO.pointlist[n+1], tetIO.pointlist[n+2]);
}

const double Mesh::get_x(int i) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    return tetIO.pointlist[n_coordinates * i + 0];
}

const double Mesh::get_y(int i) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    return tetIO.pointlist[n_coordinates * i + 1];
}

const double Mesh::get_z(int i) {
    expect(i >= 0 && i < get_n_nodes(), "Invalid index!");
    return tetIO.pointlist[n_coordinates * i + 2];
}

const double Mesh::get_node(int i, int xyz) {
    require(i >= 0 && i < get_n_nodes(), "Invalid index!");
    require(xyz >= 0 && xyz < n_coordinates, "Invalid coordinate!");
    return tetIO.pointlist[n_coordinates * i + xyz];
}

const int Mesh::get_face(int i, int node) {
    require(i >= 0 && i < get_n_faces(), "Invalid index!");
    require(node >= 0 && node < n_nodes_per_face, "Invalid node!");
    return tetIO.trifacelist[n_nodes_per_face * i + node];
}

const int Mesh::get_elem(int i, int node) {
    require(i >= 0 && i < get_n_elems(), "Invalid index!");
    require(node >= 0 && node < n_nodes_per_elem, "Invalid element!");
    return tetIO.tetrahedronlist[n_nodes_per_elem * i + node];
}

const double Mesh::get_face_centre(int i, int xyz) {
    require(i >= 0 && i < get_n_faces(), "Invalid index!");
    require(xyz >= 0 && xyz < n_coordinates, "Invalid coordinate!");

    double node1 = get_node(get_face(i, 0), xyz);
    double node2 = get_node(get_face(i, 1), xyz);
    double node3 = get_node(get_face(i, 2), xyz);
    return (node1 + node2 + node3) / 3.0;
}

const double Mesh::get_elem_centre(int i, int xyz) {
    require(i >= 0 && i < get_n_elems(), "Invalid index!");
    require(xyz >= 0 && xyz < n_coordinates, "Invalid coordinate!");

    double node1 = get_node(get_elem(i, 0), xyz);
    double node2 = get_node(get_elem(i, 1), xyz);
    double node3 = get_node(get_elem(i, 2), xyz);
    double node4 = get_node(get_elem(i, 3), xyz);
    return (node1 + node2 + node3 + node4) / 4.0;
}

const double Mesh::get_volume(const int i) {
    require(i >= 0 && i < get_n_volumes(), "Invalid index!");
    return volumes[i];
}

const int Mesh::get_nodemarker(const int i) {
    require(i >= 0 && i < get_n_nodemarkers(), "Invalid index!");
//    return tetIO.pointmarkerlist[i];
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
    return tetIO.pointlist;
}

const int* Mesh::get_faces() {
    return tetIO.trifacelist;
}

const int* Mesh::get_elems() {
    return tetIO.tetrahedronlist;
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
    tetIO.numberofpoints = N;
    tetIO.pointlist = new REAL[3 * N];
}

const void Mesh::init_faces(const int N) {
    require(N > 0, "Invalid number of faces!");
    i_faces = 0;
    tetIO.numberoftrifaces = N;
    tetIO.trifacelist = new int[3 * N];
}

const void Mesh::init_elems(const int N) {
    require(N > 0, "Invalid number of elements!");
    i_elems = 0;
    tetIO.numberoftetrahedra = N;
    tetIO.tetrahedronlist = new int[4 * N];
}

const void Mesh::init_volumes(const int N) {
    require(N > 0, "Invalid number of volumes!");
    volumes.reserve(N);
}

// =================================
// *** ADDERS: ***************
const void Mesh::add_node(const double x, const double y, const double z) {
    require(get_n_nodes() < tetIO.numberofpoints, "Allocated size of elements exceeded!");
    int i = 3 * i_nodes;
    tetIO.pointlist[i + 0] = (REAL) x;
    tetIO.pointlist[i + 1] = (REAL) y;
    tetIO.pointlist[i + 2] = (REAL) z;
    i_nodes++;
}

const void Mesh::add_face(const int f1, const int f2, const int f3) {
    require(get_n_faces() < tetIO.numberoftrifaces, "Allocated size of elements exceeded!");
    int i = 3 * i_faces;
    tetIO.trifacelist[i + 0] = f1;
    tetIO.trifacelist[i + 1] = f2;
    tetIO.trifacelist[i + 2] = f3;
    i_faces++;
}

const void Mesh::add_elem(const int e1, const int e2, const int e3, const int e4) {
    require(get_n_elems() < tetIO.numberoftetrahedra, "Allocated size of elements exceeded!");
    int i = 4 * i_elems;
    tetIO.tetrahedronlist[i + 0] = e1;
    tetIO.tetrahedronlist[i + 1] = e2;
    tetIO.tetrahedronlist[i + 2] = e3;
    tetIO.tetrahedronlist[i + 3] = e4;
    i_elems++;
}

const void Mesh::add_volume(const double V) {
    expect(get_n_volumes() < this->volumes.capacity(), "Allocated size of volumes exceeded!");
    volumes.push_back(V);
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
}

const void Mesh::copy_nodes(Mesh* mesh) {
    int N = mesh->get_n_nodes();
    for (int i = 0; i < 3 * N; ++i)
        tetIO.pointlist[i] = mesh->get_nodes()[i];
    i_nodes = N;
}

const void Mesh::copy_faces(Mesh* mesh, const int offset) {
    int N = mesh->get_n_faces();
    if (offset == 0)
        for (int i = 0; i < 3 * N; ++i)
            tetIO.trifacelist[i] = mesh->get_faces()[i];
    else
        for (int i = 0; i < 3 * N; ++i)
            tetIO.trifacelist[i] = offset + mesh->get_faces()[i];
    i_faces = N;
}

const void Mesh::copy_elems(Mesh* mesh, const int offset) {
    int N = mesh->get_n_elems();
    if (offset == 0)
        for (int i = 0; i < 4 * N; ++i)
            tetIO.tetrahedronlist[i] = mesh->get_elems()[i];
    else
        for (int i = 0; i < 4 * N; ++i)
            tetIO.tetrahedronlist[i] = offset + mesh->get_elems()[i];
    i_elems = N;
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

const void Mesh::calc_volumes() {
    int i, j, k, n1, n2, n3, n4;
    double V;
    double u[3], v[3], w[3];

    int N = get_n_elems();
    const double* nodes = get_nodes();
    const int* elems = get_elems();

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = n_nodes_per_elem * i;
        // Loop through x, y and z coordinates
        for (k = 0; k < n_coordinates; ++k) {
            n1 = n_coordinates * elems[j + 0] + k; // index of x, y or z coordinate of 1st node
            n2 = n_coordinates * elems[j + 1] + k; // ..2nd node
            n3 = n_coordinates * elems[j + 2] + k; // ..3rd node
            n4 = n_coordinates * elems[j + 3] + k; // ..4th node

            u[k] = nodes[n1] - nodes[n2];
            v[k] = nodes[n1] - nodes[n3];
            w[k] = nodes[n1] - nodes[n4];
        }
        // Volume of a rhomboedron equals to the vector product of the vectors forming its sides
        V = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]);
        // Volume of tetrahedron = 1/6 of its rhomboedron
        add_volume(fabs(V) / 6.0);
    }
}

const void Mesh::calc_volume_statistics() {
    size_t n_volumes = get_n_volumes();
    // Make a copy of the volumes vector
    vector<double> tempvec;
    tempvec.reserve(n_volumes);
    copy(volumes.begin(), volumes.end(), back_inserter(tempvec));

    sort(tempvec.begin(), tempvec.end());

    // Store min, max and average volume
    stat.Vmin = tempvec[0];
    stat.Vmax = tempvec[n_volumes - 1];
    stat.Vaverage = accumulate(tempvec.begin(), tempvec.end(), 0) / n_volumes;

    // Store the median volume
    if (n_volumes % 2 == 0)
        stat.Vmedian = (tempvec[n_volumes / 2 - 1] + tempvec[n_volumes / 2]) / 2;
    else
        stat.Vmedian = tempvec[n_volumes / 2];
}

// Function to perform tetgen calculation on input and output data
const void Mesh::recalc(const string cmd) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, &tetIO, &tetIO);

    i_nodes = tetIO.numberofpoints;
    i_faces = tetIO.numberoftrifaces;
    i_elems = tetIO.numberoftetrahedra;
}

// Write tetgen mesh into files with its internal functions
// Files are written into working directory with their default names
const void Mesh::output() {
    const string cmd = "Q";
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, &tetIO, NULL);
}

// Function to output mesh in .vtk format
const void Mesh::write_vtk(const string file_name, const int n_nodes, const int n_cells,
        const int n_markers, const REAL* nodes, const int* cells, const vector<int>* markers,
        const int celltype, const int n_nodes_in_cell) {

    int i, j;
    char file_name_char[1024];
    strcpy(file_name_char, file_name.c_str());

    FILE *out_file;
    out_file = fopen(file_name_char, "w");
    require(out_file != (FILE *) NULL, "Can't open a file " + file_name);

//    fprintf(outfile, "# vtk DataFile Version 2.0\n");
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

    // Output the cells (tetrahedra or triangles)
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

    // Output the types of cells, 10=tetrahedron, 5=triangle
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

// Function to output nodes in .xyz format
const void Mesh::write_xyz(const string file_name, const int n_nodes, const int n_markers, const REAL* nodes, const vector<int>* markers) {
    ofstream out_file(file_name);
    require(out_file.is_open(), "Can't open a file " + file_name);

    out_file << n_nodes << "\n";
    out_file << "Data of the nanotip: id x y z marker\n";

    for (int i = 0; i < n_nodes; ++i) {
        out_file << i << " ";
        out_file << nodes[3*i+0] << " ";
        out_file << nodes[3*i+1] << " ";
        out_file << nodes[3*i+2] << " ";

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

    int n_nodes = get_n_nodes();
    int n_markers = get_n_nodemarkers();
    const REAL* nodes = get_nodes();           // pointer to nodes
    vector<int>* markers = &nodemarkers; // pointer to node markers

    if (file_type == "vtk") {
        const int celltype = 1;     // 1-vertex, 5-triangle, 10-tetrahedron

        int node_indx[n_nodes];
        for(int i = 0; i < n_nodes; ++i)
            node_indx[i] = i;

        write_vtk(file_name, n_nodes, n_nodes, n_markers, nodes, node_indx, markers, celltype, 1);

    } else
        write_xyz(file_name, n_nodes, n_markers, nodes, markers);

}

// Function to extract file type from file name
const string Mesh::get_file_type(const string file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    return file_name.substr(start, end);
}

// =================================
} /* namespace femocs */
