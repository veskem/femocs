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

void Mesher::add_vacuum(const string cmd, Femocs::SimuCell* cell) {
    tetgenio* tetIO = &tetgenOut;

    int i, j, N, M;

    N = tetIO->numberofpoints;
    M = N + 4;
    vector<double> pointstore;
    pointstore.reserve(3*N);
    for (i = 0; i < 3*N; ++i)
        pointstore.push_back(tetIO->pointlist[i]);

    vector<int> pointmarkerstore;
    pointmarkerstore.reserve(N);
    for (i = 0; i < N; ++i)
        pointmarkerstore.push_back(tetIO->pointmarkerlist[i]);

    tetIO->numberofpoints = M;
    tetIO->pointlist = new REAL[3 * M];
    tetIO->pointmarkerlist = new int[M];

    for (i = 0; i < 3*N; ++i)
        tetIO->pointlist[i] = pointstore[i];

    for (i = 0; i < N; ++i)
        tetIO->pointmarkerlist[i] = pointmarkerstore[i];

    j = 3*N;
    add_point(cell->xmin, cell->ymin, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmin, cell->ymax, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmax, cell->ymin, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmax, cell->ymax, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;

    calculate(cmd, tetIO, tetIO);
}
// Function to generate mesh between surface atoms
void Mesher::generate_mesh(const Femocs::SimuCell* cell, shared_ptr<Surface> bulk,
        shared_ptr<Surface> surf, Vacuum* vacuum, string cmd) {

//    generate_surface_mesh(surf);
    generate_bulk_mesh(cell, bulk);
//    generate_volume_mesh(bulk, vacuum);
    calculate(cmd, &tetgenOut, &tetgenOut);
}

void Mesher::unite3(const tetgenio* tetIO_bulk, const tetgenio* tetIO_surf, const Femocs::SimuCell* cell, const string cmd) {
    tetgenio* tetIO_new = &tetgenOut;

    int i, j, k;

    int N_elems = tetIO_bulk->numberoftetrahedra;
    int N_points = tetIO_bulk->numberofpoints;
    int N_surf = tetIO_surf->numberoftrifaces;
    int N_bulk = tetIO_bulk->numberoftrifaces;
    vector<bool> isQuality(N_surf);

    // Loop through the faces
    for (i = 0; i < N_surf; ++i) {
        isQuality[i] = tetIO_surf->trifacemarkerlist[i] == cell->type_none;
    }

    int M = N_bulk + accumulate(isQuality.begin(), isQuality.end(), 0);

    // Copy the points of bulk without modification
    tetIO_new->numberofpoints = N_points;
    tetIO_new->pointlist = new REAL[3*N_points];
    for (i = 0; i < 3*N_points; ++i)
        tetIO_new->pointlist[i] = tetIO_bulk->pointlist[i];

    // Copy the elements of bulk without modification
    tetIO_new->numberoftetrahedra = N_elems;
    tetIO_new->tetrahedronlist = new int[4*N_elems];
    for (i = 0; i < 4*N_elems; ++i)
        tetIO_new->tetrahedronlist[i] = tetIO_bulk->tetrahedronlist[i];

    // Copy the faces of bulk without modification
    tetIO_new->numberoftrifaces = M;
    tetIO_new->trifacelist = new int[3*M];
    for (i = 0; i < 3*N_bulk; ++i)
        tetIO_new->trifacelist[i] = tetIO_bulk->trifacelist[i];

    // Copy only the surface faces that are not on the edge of simulation cell
    j = 3*N_bulk;
    for (i = 0; i < N_surf; ++i) {
        if ( isQuality[i] ) {
            k = 3*i;
            tetIO_new->trifacelist[j+0] = tetIO_surf->trifacelist[k+0];
            tetIO_new->trifacelist[j+1] = tetIO_surf->trifacelist[k+1];
            tetIO_new->trifacelist[j+2] = tetIO_surf->trifacelist[k+2];
            j += 3;
        }
    }

    calculate(cmd, tetIO_new, tetIO_new);
}

void Mesher::unite2(const string cmd, const Mesher* mesh2, const Femocs::SimuCell* cell) {
    tetgenio* tetIO = &tetgenOut;
    REAL* node = mesh2->tetgenOut.pointlist;      // pointer to face mesh nodes
    int* face = mesh2->tetgenOut.trifacelist;     // pointer to face mesh faces

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };

    int i, j, k, l1, l2, l3;
    bool on_simucell_edge;

    int N = tetIO->numberoftrifaces;
    int N2 = mesh2->tetgenOut.numberoftrifaces;
    vector<bool> isQuality(N2);

    // Loop through the faces
    for (i = 0; i < N2; ++i) {
        j = 3 * i;
        on_simucell_edge = 0;
        // Loop through x, y and z coordinates
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k; // ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            // Get whether the face is on the boundary of simulation cell
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2*k+0]);
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2*k+1]);
        }
        // Keep only the faces with appropriate quality
        isQuality[i] = 1;//!on_simucell_edge;
    }

    int M = N + accumulate(isQuality.begin(), isQuality.end(), 0);

    // Save current mesher faces with their markers into temporary array
    vector<int> facelistsave;
    facelistsave.reserve(3*N);
    for (i = 0; i < 3*N; ++i)
        facelistsave.push_back(tetIO->trifacelist[i]);

    // Make new facelist and first fill it with old values
    tetIO->numberoftrifaces = M;
    tetIO->trifacelist = new int[3*M];
    for (i = 0; i < 3*N; ++i)
        tetIO->trifacelist[i] = facelistsave[i];

    // ... and then add faces from parameter mesh that are not on the edge of simulation cell
    j = 3*N;
    for (i = 0; i < N2; ++i) {
        if ( isQuality[i] ) {
            k = 3*i;
            tetIO->trifacelist[j+0] = face[k+0];
            tetIO->trifacelist[j+1] = face[k+1];
            tetIO->trifacelist[j+2] = face[k+2];
            j += 3;
        }
    }

    calculate(cmd, tetIO, tetIO);
}

void Mesher::unite(const string cmd, const Mesher* mesh2) {
    tetgenio* tetIO = &tetgenOut;
    int N = tetIO->numberoftrifaces;
    int N2 = mesh2->tetgenOut.numberoftrifaces;
    int M = N + N2;

    int i, j;
    vector<int> facelistsave;
    facelistsave.reserve(3*N);
    for (i = 0; i < 3*N; ++i)
        facelistsave.push_back(tetIO->trifacelist[i]);

    tetIO->numberoftrifaces = M;
    tetIO->trifacelist = new int[3*M];

    for (i = 0; i < 3*N; ++i)
        tetIO->trifacelist[i] = facelistsave[i];

    for (i = 0; i < 3*N2; ++i) {
        j = 3*N + i;
        tetIO->trifacelist[j] = mesh2->tetgenOut.trifacelist[i];
    }
    calculate(cmd, tetIO, tetIO);
}

// Function to generate volumetric mesh between surface atoms and vacuum
void Mesher::generate_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum) {
    int N_bulk = bulk->getN();      // number of atoms in bulk material
    int N_vacuum = vacuum->getN();  // number of atoms in vacuum
    int M = N_bulk + N_vacuum;      // total number of nodes

    tetgenio* tetIO = &tetgenOut;
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

// Function to generate mesh between surface atoms
void Mesher::generate_bulk_mesh(const Femocs::SimuCell* cell, shared_ptr<Surface> bulk) {
    int i, j;
    int N = bulk->getN();
    int M = N+4;

    tetgenio* tetIO = &tetgenOut;
    tetIO->numberofpoints = M;
    tetIO->pointlist = new REAL[3 * M];
    tetIO->pointmarkerlist = new int[M];

    j = 0;
    for (i = 0; i < N; ++i) {
        add_point(bulk->getX(i), bulk->getY(i), bulk->getZ(i), bulk->getType(i), j, tetIO);
        j += 3;
    }
    add_point(cell->xmin, cell->ymin, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmin, cell->ymax, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmax, cell->ymin, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
    add_point(cell->xmax, cell->ymax, cell->zmaxbox, cell->type_vacuum, j, tetIO);
    j += 3;
}

// Function to generate mesh between surface atoms
void Mesher::generate_surface_mesh(shared_ptr<Surface> surf) {
    int i, j;
    int N = surf->getN();

    tetgenio* tetIO = &tetgenOut;
    tetIO->numberofpoints = N;
    tetIO->pointlist = new REAL[3 * N];
    tetIO->pointmarkerlist = new int[N];

    for (i = 0; i < N; ++i) {
        j = 3 * i;
        add_point(surf->getX(i), surf->getY(i), surf->getZ(i), surf->getType(i), j, tetIO);
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
    const int markers[6] = { cell->type_xmin, cell->type_xmax, cell->type_ymin, cell->type_ymax,
            cell->type_zmin, cell->type_zmax };

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        emarker[i] = cell->type_none; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * elem[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * elem[j + 1] + k; // ..2nd node
            l3 = ncoords * elem[j + 2] + k; // ..3rd node
            l4 = ncoords * elem[j + 3] + k; // ..4th node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // element on negative direction edge
            if ( on_face(node[l1], node[l2], node[l3], xyz[m1])
                    || on_face(node[l1], node[l2], node[l4], xyz[m1])
                    || on_face(node[l1], node[l3], node[l4], xyz[m1])
                    || on_face(node[l2], node[l3], node[l4], xyz[m1]) )
                emarker[i] = markers[m1];

            // element on positive direction edge
            if ( on_face(node[l1], node[l2], node[l3], xyz[m2])
                    || on_face(node[l1], node[l2], node[l4], xyz[m2])
                    || on_face(node[l1], node[l3], node[l4], xyz[m2])
                    || on_face(node[l2], node[l3], node[l4], xyz[m2]) )
                emarker[i] = markers[m2];
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
    const int markers[6] = { cell->type_xmin, cell->type_xmax, cell->type_ymin, cell->type_ymax,
            cell->type_zmin, cell->type_zmax };

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 3;  // nr of nodes per face

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        fmarker[i] = cell->type_none; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * face[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * face[j + 1] + k; // ..2nd node
            l3 = ncoords * face[j + 2] + k; // ..3rd node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // outer face in negative direction
            if (on_face(node[l1], node[l2], node[l3], xyz[m1]))
                fmarker[i] = markers[m1];

            // outer face in positive direction
            if (on_face(node[l1], node[l2], node[l3], xyz[m2]))
                fmarker[i] = markers[m2];
        }
    }
}

// Determine whether the nodes of two faces are overlapping
bool Mesher::on_face(const double f1n1, const double f1n2, const double f1n3, const double f2) {
    double eps = 0.1;
    bool dif1 = (fabs(f1n1 - f2) < eps);
    bool dif2 = (fabs(f1n2 - f2) < eps);
    bool dif3 = (fabs(f1n3 - f2) < eps);
    return dif1 && dif2 && dif3;
}

// Function to clean the mesh from hull elements
void Mesher::clean_mesh(string cmd, double rmax, const Femocs::SimuCell* cell) {
    clean_elements(&tetgenOut, rmax);
    clean_faces2(&tetgenOut, rmax, cell);
//    clean_faces(&tetgenOut, rmax);
    clean_edges(&tetgenOut, rmax);
    calculate(cmd, &tetgenOut, &tetgenOut);
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
void Mesher::clean_faces2(tetgenio* tetIO, double rmax, const Femocs::SimuCell* cell) {
    REAL* node = tetIO->pointlist;      // pointer to nodes
    int* face = tetIO->trifacelist;     // pointer to triangular faces

    double dx, r12, r13, r23;
    int i, j, k, l1, l2, l3;
    bool on_simucell_edge;
    int N = tetIO->numberoftrifaces;
    vector<bool> isQuality(N);

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r23 = 0;
        on_simucell_edge = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k;   // ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            dx = node[l1] - node[l2];
            r12 += dx * dx; // length^2 of face's 1st edge
            dx = node[l1] - node[l3];
            r13 += dx * dx; // ..2nd edge
            dx = node[l2] - node[l3];
            r23 += dx * dx; // ..3rd edge

            // Detect whether the face is on the edge of simulation box
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2*k+0]);
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2*k+1]);
        }
        // Keep only the faces with appropriate quality
        isQuality[i] = (r12 < rmax) & (r13 < rmax) & (r23 < rmax) & (!on_simucell_edge);
    }
    // Insert faces with suitable quality
    update_list(tetIO->trifacelist, &(tetIO->numberoftrifaces), isQuality, 3);
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

// Public function to perform tetgen calculation cycle
void Mesher::calc(const string cmd) {
    calculate(cmd, &tetgenOut, &tetgenOut);
}

} /* namespace femocs */
