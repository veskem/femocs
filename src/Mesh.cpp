/*
 * Mesh.cpp
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#include "Mesh.h"

#include <stdio.h>
#include <cstring>

using namespace std;
namespace femocs {

Mesh::Mesh() {
    inodes = 0;
    ielems = 0;
    ifaces = 0;
}

double* Mesh::getNodes() {
    return tetIO.pointlist;
}

int* Mesh::getFaces() {
    return tetIO.trifacelist;
}

int* Mesh::getElems() {
    return tetIO.tetrahedronlist;
}

const int Mesh::getNodemarker(const int i) {
    return nodemarkers[i];
}

const int Mesh::getFacemarker(const int i) {
    return facemarkers[i];
}

const int Mesh::getElemmarker(const int i) {
    return elemmarkers[i];
}

const int Mesh::getNnodes() {
    return tetIO.numberofpoints;
}

const int Mesh::getNelems() {
    return tetIO.numberoftetrahedra;
}

const int Mesh::getNfaces() {
    return tetIO.numberoftrifaces;
}

const int Mesh::getNnodemarkers() {
    return nodemarkers.size();
}

const int Mesh::getNfacemarkers() {
    return facemarkers.size();
}

const int Mesh::getNelemmarkers() {
    return elemmarkers.size();
}

void Mesh::init_facemarkers(const int N) {
    facemarkers.reserve(N);
}

void Mesh::init_elemmarkers(const int N) {
    elemmarkers.reserve(N);
}

void Mesh::init_nodemarkers(const int N) {
    nodemarkers.reserve(N);
}

void Mesh::init_nodes(const int N) {
    inodes = 0;
    tetIO.numberofpoints = N;
    tetIO.pointlist = new REAL[3 * N];
}

void Mesh::init_faces(const int N) {
    ifaces = 0;
    tetIO.numberoftrifaces = N;
    tetIO.trifacelist = new int[3 * N];
}

void Mesh::init_elems(const int N) {
    ielems = 0;
    tetIO.numberoftetrahedra = N;
    tetIO.tetrahedronlist = new int[4 * N];
}

void Mesh::add_nodemarker(const int m) {
    nodemarkers.push_back(m);
}

void Mesh::add_facemarker(const int m) {
    facemarkers.push_back(m);
}

void Mesh::add_elemmarker(const int m) {
    elemmarkers.push_back(m);
}

void Mesh::add_elem(const int e1, const int e2, const int e3, const int e4) {
    int i = 4 * ielems;
    tetIO.tetrahedronlist[i + 0] = e1;
    tetIO.tetrahedronlist[i + 1] = e2;
    tetIO.tetrahedronlist[i + 2] = e3;
    tetIO.tetrahedronlist[i + 4] = e4;
    ielems++;
}

void Mesh::add_face(const int f1, const int f2, const int f3) {
    int i = 3 * ifaces;
    tetIO.trifacelist[i + 0] = f1;
    tetIO.trifacelist[i + 1] = f2;
    tetIO.trifacelist[i + 2] = f3;
    ifaces++;
}

void Mesh::add_node(const double x, const double y, const double z) {
    int i = 3 * inodes;
    tetIO.pointlist[i + 0] = (REAL) x;
    tetIO.pointlist[i + 1] = (REAL) y;
    tetIO.pointlist[i + 2] = (REAL) z;
    inodes++;
}

void Mesh::copy_nodes(shared_ptr<Mesh> mesh) {
    int N = mesh->getNnodes();
    for (int i = 0; i < 3 * N; ++i)
        tetIO.pointlist[i] = mesh->getNodes()[i];
    inodes = N;
}

void Mesh::copy_faces(shared_ptr<Mesh> mesh) {
    int N = mesh->getNfaces();
    for (int i = 0; i < 3 * N; ++i)
        tetIO.trifacelist[i] = mesh->getFaces()[i];
    ifaces = N;
}

void Mesh::copy_elems(shared_ptr<Mesh> mesh) {
    int N = mesh->getNelems();
    for (int i = 0; i < 4 * N; ++i)
        tetIO.tetrahedronlist[i] = mesh->getElems()[i];
    ielems = N;
}

// Function to perform tetgen calculation on input and output data
void Mesh::calculate(string cmd) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, &tetIO, &tetIO);
}

// Function to output mesh in .vtk format
void Mesh::write_vtk(const string file_name, const int nnodes, const int ncells, const int nmarkers,
        const REAL* nodes, const int* cells, const vector<int>* markers, const int celltype,
        const int nnodes_per_cell) {
    int i, j;
    char file_name_char[1024];
    strcpy(file_name_char, file_name.c_str());

    FILE *outfile;
    outfile = fopen(file_name_char, "w");
    if (outfile == (FILE *) NULL) {
        printf("File I/O Error:  Cannot create file %s.\n", file_name_char);
        return;
    }

    fprintf(outfile, "# vtk DataFile Version 2.0\n");
    fprintf(outfile, "Unstructured Grid\n");
    fprintf(outfile, "ASCII\n"); // another option is BINARY
    fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");

//    fprintf(outfile, "# vtk DataFile Version 3.0\n");
//    fprintf(outfile, "# This file was generated for test purposes\n");
//    fprintf(outfile, "ASCII\n"); // another option is BINARY
//    fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n\n");

    // Output the nodes
    if (nnodes > 0) {
        fprintf(outfile, "POINTS %d double\n", nnodes);
        for (i = 0; i < 3 * nnodes; i += 3)
            fprintf(outfile, "%.8g %.8g %.8g\n", nodes[i + 0], nodes[i + 1], nodes[i + 2]);
    }

    // Output the cells (tetrahedra or triangles)
    if (ncells > 0) {
        fprintf(outfile, "CELLS %d %d\n", ncells, ncells * (nnodes_per_cell + 1));
        for (i = 0; i < nnodes_per_cell * ncells; i += nnodes_per_cell) {
            fprintf(outfile, "%d ", nnodes_per_cell);
            for (j = 0; j < nnodes_per_cell; ++j)
                fprintf(outfile, "%d ", cells[i + j]);
            fprintf(outfile, "\n");
        }
    }

    // Output the types of cells, 10=tetrahedron, 5=triangle
    if (ncells > 0) {
        fprintf(outfile, "CELL_TYPES %d\n", ncells);
        for (i = 0; i < ncells; ++i)
            fprintf(outfile, "%d\n", celltype);
    }

    // Output cell attributes
    if (nmarkers > 0) {
        fprintf(outfile, "CELL_DATA %d\n", nmarkers);
        fprintf(outfile, "SCALARS Cell_markers int\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < nmarkers; ++i)
            fprintf(outfile, "%d\n", facemarkers[i]);
    }

    fclose(outfile);
}

// Function to output faces in .vtk format
void Mesh::write_faces(const string file_name) {
    const int celltype = 5;
    const int nnodes_in_cell = 3;

    int nnodes = tetIO.numberofpoints;
    int nfaces = tetIO.numberoftrifaces;
    int nmarkers = getNfacemarkers();
    REAL* nodes = tetIO.pointlist;          // pointer to nodes
    int* faces = tetIO.trifacelist;         // pointer to face nodes
    vector<int>* markers = &facemarkers;    // pointer to face markers

    write_vtk(file_name, nnodes, nfaces, nmarkers, nodes, faces, markers, celltype, nnodes_in_cell);
}

// Function to output faces in .vtk format
void Mesh::write_elems(const string file_name) {
    const int celltype = 10; // 5-triangle, 10-tetrahedron
    const int nnodes_in_cell = 4;

    int nnodes = tetIO.numberofpoints;
    int nelems = tetIO.numberoftetrahedra;
    int nmarkers = getNelemmarkers();
    REAL* nodes = tetIO.pointlist;          // pointer to nodes
    int* elems = tetIO.tetrahedronlist;     // pointer to faces nodes
    vector<int>* markers = &elemmarkers;    // pointer to element markers

    write_vtk(file_name, nnodes, nelems, nmarkers, nodes, elems, markers, celltype, nnodes_in_cell);
}

} /* namespace femocs */
