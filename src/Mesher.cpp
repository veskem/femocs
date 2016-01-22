/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

#include <stdio.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <sstream>

#include "Vacuum.h"

using namespace std;
namespace femocs {

Mesher::Mesher(string mesher) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
}

// Function to generate mesh between surface atoms
void Mesher::generate_mesh(const Femocs::SimuCell* cell, shared_ptr<Surface> bulk,
        shared_ptr<Surface> surf, Vacuum* vacuum, string cmd) {
//    generate_surface_mesh(surf, &tetgenIn);
    generate_volume_mesh(bulk, vacuum, &tetgenIn);
    calculate(cmd, &tetgenIn, &tetgenOut);
}

// Function to generate volumetric mesh between surface atoms and vacuum
void Mesher::generate_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum, tetgenio* tetIO) {
    int N_bulk = bulk->getN();      // number of atoms in bulk material
    int N_vacuum = vacuum->getN();  // number of atoms in vacuum
    int M = N_bulk + N_vacuum;      // total number of nodes

    tetIO->numberofpoints = M;
    tetIO->pointlist = new REAL[3 * M];
    tetIO->pointmarkerlist = new int[M];

    int i;
    int j = 0;
    // Add bulk atoms
    for (i = 0; i < N_bulk; ++i) {
        add_point(bulk->getX(i), bulk->getY(i), bulk->getZ(i), bulk->getType(i), j, tetIO);
        j += 3;
    }
    // Add vacuum atoms
    for (i = 0; i < N_vacuum; ++i) {
        add_point(vacuum->getX(i), vacuum->getY(i), vacuum->getZ(i), vacuum->getType(i), j, tetIO);
        j += 3;
    }
}

// Function to add point with coordinates and maker to tetgen list
void Mesher::add_point(const double x, const double y, const double z, const int pmarker,
        const int i, tetgenio* tetIO) {
    tetIO->pointmarkerlist[i / 3] = pmarker;
    tetIO->pointlist[i + 0] = (REAL) x;
    tetIO->pointlist[i + 1] = (REAL) y;
    tetIO->pointlist[i + 2] = (REAL) z;
}

// Function to generate mesh between surface atoms
void Mesher::generate_surface_mesh(shared_ptr<Surface> surf, tetgenio* tetIO) {
    int i, j;
    int N = surf->getN();

    tetIO->numberofpoints = N;
    tetIO->pointlist = new REAL[3 * N];
    tetIO->pointmarkerlist = new int[N];

    for (i = 0; i < N; ++i) {
        j = 3 * i;
        tetIO->pointlist[j + 0] = (REAL) surf->getX(i);
        tetIO->pointlist[j + 1] = (REAL) surf->getY(i);
        tetIO->pointlist[j + 2] = (REAL) surf->getZ(i);
        tetIO->pointmarkerlist[i] = surf->getType(i);
    }
}

// Function to clean the mesh from hull elements
void Mesher::clean_mesh(string cmd, double rmax) {
    clean_elements(&tetgenOut, rmax);
    clean_faces(&tetgenOut, rmax);
    clean_edges(&tetgenOut, rmax);
    calculate(cmd, &tetgenOut, &tetgenOut);
}

// Function to mark the elements and faces of mesh
void Mesher::mark_mesh(const Femocs::SimuCell* cell, string cmd) {
    mark_faces(cell, &tetgenOut);
    mark_elems(cell, &tetgenOut);
    calculate(cmd, &tetgenOut, &tetgenOut);
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_elems(const Femocs::SimuCell* cell, tetgenio* tetIO) {
    int i, j, k, l1, l2, l3, l4, m1, m2;
    int N = tetIO->numberoftetrahedra;

    REAL* node = tetIO->pointlist;       // pointer to nodes
    int* elem = tetIO->tetrahedronlist;  // pointer to tetrahedrons
    tetIO->tetrahedronattributelist = new double[N];
    double* emarker = tetIO->tetrahedronattributelist; // pointer to tetrahedron attributes

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };
    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        emarker[i] = ncoords; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * elem[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * elem[j + 1] + k; // ..2nd node
            l3 = ncoords * elem[j + 2] + k; // ..3rd node
            l4 = ncoords * elem[j + 3] + k; // ..4th node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // element on negative direction edge
            if (on_face(node[l1], node[l2], node[l3], xyz[m1], xyz[m1], xyz[m1])
                    || on_face(node[l1], node[l2], node[l4], xyz[m1], xyz[m1], xyz[m1])
                    || on_face(node[l1], node[l3], node[l4], xyz[m1], xyz[m1], xyz[m1])
                    || on_face(node[l2], node[l3], node[l4], xyz[m1], xyz[m1], xyz[m1]))
                emarker[i] = ncoords - k - 1;

            // element on positive direction edge
            if (on_face(node[l1], node[l2], node[l3], xyz[m2], xyz[m2], xyz[m2])
                    || on_face(node[l1], node[l2], node[l4], xyz[m2], xyz[m2], xyz[m2])
                    || on_face(node[l1], node[l3], node[l4], xyz[m2], xyz[m2], xyz[m2])
                    || on_face(node[l2], node[l3], node[l4], xyz[m2], xyz[m2], xyz[m2]))
                emarker[i] = ncoords + k + 1;
        }
    }
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_faces(const Femocs::SimuCell* cell, tetgenio* tetIO) {
    int i, j, k, l1, l2, l3, m1, m2;
    int N = tetIO->numberoftrifaces;

    REAL* node = tetIO->pointlist;	// pointer to nodes
    int* face = tetIO->trifacelist;	// pointer to tetrahedron faces
    tetIO->trifacemarkerlist = new int[N];
    int* fmarker = tetIO->trifacemarkerlist; // pointer to face markers

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };
    const int ncoords = 3; // nr of coordinates
    const int nnodes = 3;  // nr of nodes per face

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        fmarker[i] = ncoords; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * face[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * face[j + 1] + k; // ..2nd node
            l3 = ncoords * face[j + 2] + k; // ..3rd node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // outer face in negative direction
            if (on_face(node[l1], node[l2], node[l3], xyz[m1], xyz[m1], xyz[m1]))
                fmarker[i] = ncoords - k - 1;

            // outer face in positive direction
            if (on_face(node[l1], node[l2], node[l3], xyz[m2], xyz[m2], xyz[m2]))
                fmarker[i] = ncoords + k + 1;
        }
    }
}

// Determine whether the nodes of two faces are overlapping
bool Mesher::on_face(const double f1n1, const double f1n2, const double f1n3, const double f2n1,
        const double f2n2, const double f2n3) {
    double eps = 0.1;
    bool dif1 = (fabs(f1n1 - f2n1) < eps);
    bool dif2 = (fabs(f1n2 - f2n2) < eps);
    bool dif3 = (fabs(f1n3 - f2n3) < eps);
    return dif1 && dif2 && dif3;
}

// Function to remove too big tetrahedra from the mesh
void Mesher::clean_elements(tetgenio* tetIO, double rmax) {
    REAL* node = tetIO->pointlist;		// pointer to nodes
    int* elem = tetIO->tetrahedronlist;	// pointer to tetrahedron corners

    double dx, r12, r13, r14, r23, r24, r34;
    int i, j, k, l1, l2, l3, l4;

    int N = tetIO->numberoftetrahedra;
    vector<bool> isQuality(N);

    // Loop through the tetrahedra
    for (i = 0; i < N; ++i) {
        j = 4 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r14 = r23 = r24 = r34 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * elem[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * elem[j + 1] + k;	// ..2nd node
            l3 = 3 * elem[j + 2] + k; // ..3rd node
            l4 = 3 * elem[j + 3] + k; // ..4th node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of tetrahedron 1st edge
            dx = node[l1] - node[l3];
            r13 += dx * dx; // ..2nd edge
            dx = node[l1] - node[l4];
            r14 += dx * dx; // ..3rd edge
            dx = node[l2] - node[l3];
            r23 += dx * dx; // ..4th edge
            dx = node[l2] - node[l4];
            r24 += dx * dx; // ..5th edge
            dx = node[l3] - node[l4];
            r34 += dx * dx; // ..6th edge
        }
        // Keep only the elements with appropriate quality
        isQuality[i] = (r12 < rmax) & (r13 < rmax) & (r14 < rmax) & (r23 < rmax) & (r24 < rmax)
                & (r34 < rmax);
    }
    // Insert tetrahedra with suitable quality
    update_list(tetIO->tetrahedronlist, &(tetIO->numberoftetrahedra), isQuality, 4);
}

// Function to remove too big faces from the mesh
void Mesher::clean_faces(tetgenio* tetIO, double rmax) {
    REAL* node = tetIO->pointlist;		// pointer to nodes
    int* face = tetIO->trifacelist;		// pointer to triangular faces

    double dx, r12, r13, r23;
    int i, j, k, l1, l2, l3;

    int N = tetIO->numberoftrifaces;
    vector<bool> isQuality(N);

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r23 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k;	// ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of face's 1st edge
            dx = node[l1] - node[l3];
            r13 += dx * dx; // ..2nd edge
            dx = node[l2] - node[l3];
            r23 += dx * dx; // ..3rd edge
        }
        // Keep only the faces with appropriate quality
        isQuality[i] = (r12 < rmax) & (r13 < rmax) & (r23 < rmax);
    }
    // Insert faces with suitable quality
    update_list(tetIO->trifacelist, &(tetIO->numberoftrifaces), isQuality, 3);
}

// Function to remove too big edges from the mesh
void Mesher::clean_edges(tetgenio* tetIO, double rmax) {
    REAL* node = tetIO->pointlist;	// pointer to nodes
    int* edge = tetIO->edgelist;	// pointer to edges

    double dx, r12;
    int i, j, k, l1, l2;

    int N = tetIO->numberofedges;
    vector<bool> isQuality(N);

    // Loop through the edges
    for (i = 0; i < N; ++i) {
        j = 2 * i;
        // Loop through x, y and z coordinates
        r12 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 2 * edge[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 2 * edge[j + 1] + k;	// ..2nd node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of edge
        }
        // Keep only the edges with appropriate quality
        isQuality[i] = (r12 < rmax);
    }
    // Insert edges with suitable quality
    update_list(tetIO->edgelist, &(tetIO->numberofedges), isQuality, 2);
}

// Function to remove the objects from element, face or edge list
void Mesher::update_list(int* list, int* list_size, vector<bool> is_quality, int M) {
    int i, j, k;
    int N = is_quality.size();
    // new list size = number of objects with appropriate quality
    *list_size = accumulate(is_quality.begin(), is_quality.end(), 0);
    j = 0;
    // loop through the old array and move quality objects to the left
    for (i = 0; i < N; ++i) {
        if (is_quality[i]) {
            ++j;
            for (k = 0; k < M; ++k)
                list[M * j + k] = list[M * i + k];
        }
    }
}

// Function to output mesh in .vtk format
void Mesher::write_vtk(const string file_name, const int nnodes, const int ncells,
        const int nmarkers, const REAL* nodes, const int* cells, const char* markerstr,
        const int celltype, const int nnodes_per_cell) {
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

    // Output the nodes
    if (nnodes > 0) {
        fprintf(outfile, "POINTS %d double\n", nnodes);
        for (i = 0; i < 3 * nnodes; i += 3)
            fprintf(outfile, "%.8g %.8g %.8g\n", nodes[i + 0], nodes[i + 1], nodes[i + 2]);
        fprintf(outfile, "\n");
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
        fprintf(outfile, "\n");
    }

    // Output the types of cells, 10=tetrahedron, 5=triangle
    if (ncells > 0) {
        fprintf(outfile, "CELL_TYPES %d\n", ncells);
        for (i = 0; i < ncells; ++i)
            fprintf(outfile, "%d\n", celltype);
        fprintf(outfile, "\n");
    }

    // Output cell attributes
    if (nmarkers > 0) {
        fprintf(outfile, "CELL_DATA %d\n", nmarkers);
        fprintf(outfile, "SCALARS Cell_markers int\n");
        fprintf(outfile, "LOOKUP_TABLE default\n");
        fprintf(outfile, "%s", markerstr);
    }

    fclose(outfile);
}

// Function to output faces in .vtk format
void Mesher::write_faces(const string file_name) {
    tetgenio* tetIO;
    tetIO = &tetgenOut;
    ostringstream markerstr;

    const int celltype = 5;
    const int nnodes_per_cell = 3;

    int nnodes = tetIO->numberofpoints;
    int nfaces = tetIO->numberoftrifaces;
    int nmarkers = tetIO->numberoftrifaces;
    REAL* nodes = tetIO->pointlist;              // pointer to nodes
    int* faces = tetIO->trifacelist;             // pointer to face nodes
    int* markers = tetIO->trifacemarkerlist;    // pointer to face markers

    // Convert cell attributes to string stream
    for (int i = 0; i < nmarkers; ++i)
        markerstr << markers[i] << "\n";

    write_vtk(file_name, nnodes, nfaces, nmarkers, nodes, faces, markerstr.str().c_str(), celltype,
            nnodes_per_cell);
}

// Function to output faces in .vtk format
void Mesher::write_elems(const string file_name) {
    tetgenio* tetIO;
    tetIO = &tetgenOut;
    ostringstream markerstr;

    const int celltype = 10; // 5-triangle, 10-tetrahedron
    const int nnodes_per_cell = 4;

    int nnodes = tetIO->numberofpoints;
    int nelems = tetIO->numberoftetrahedra;
    int nmarkers = tetIO->numberoftetrahedra;
    REAL* nodes = tetIO->pointlist;              // pointer to nodes
    int* elems = tetIO->tetrahedronlist;             // pointer to faces nodes
    double* markers = tetIO->tetrahedronattributelist;    // pointer to faces markers

    // Convert cell attributes to string stream
    for (int i = 0; i < nmarkers; ++i)
        markerstr << (int) markers[i] << "\n";

    write_vtk(file_name, nnodes, nelems, nmarkers, nodes, elems, markerstr.str().c_str(), celltype,
            nnodes_per_cell);
}

// Function to perform standalone mesh refinement operations
void Mesher::refine(const string cmd) {
    calculate(cmd, &tetgenOut, &tetgenOut);
}

// Function to perform tetgen calculation on input and output data
void Mesher::calculate(string cmd, tetgenio* t1, tetgenio* t2) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, t1, t2);
}

} /* namespace femocs */
