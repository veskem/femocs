/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

#include <iostream>
#include <numeric>
#include <vector>

using namespace std;
namespace femocs {

Mesher::Mesher(string mesher) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
}

// Function to generate mesh between surface atoms
void Mesher::generate_mesh(shared_ptr<Surface> surf, string cmd) {
    int i, j;
    int N = surf->getN();

    tetgenIn.numberofpoints = N;
    tetgenIn.pointlist = new REAL[3 * N];
    tetgenIn.pointmarkerlist = new int[N];

    for (i = 0; i < N; ++i) {
        j = 3 * i;
        tetgenIn.pointlist[j + 0] = (REAL) surf->getX(i);
        tetgenIn.pointlist[j + 1] = (REAL) surf->getY(i);
        tetgenIn.pointlist[j + 2] = (REAL) surf->getZ(i);
        tetgenIn.pointmarkerlist[i] = surf->getType(i);
    }

    calculate(cmd, &tetgenIn, &tetgenOut);
}

// Function to clean the mesh from hull elements
void Mesher::clean_mesh(string cmd, double rmax) {
    clean_elements(&tetgenOut, rmax);
    clean_faces(&tetgenOut, rmax);
    clean_edges(&tetgenOut, rmax);

    calculate(cmd, &tetgenOut, &tetgenOut);
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_mesh(AtomReader::Data data, string cmd) {
    int i, j, k, l1, l2, l3;
    int N = tetgenOut.numberoftrifaces;

    REAL* node = tetgenOut.pointlist;	// pointer to nodes
    int* face = tetgenOut.trifacelist;	// pointer to tetrahedron faces
    int* fmarker = tetgenOut.trifacemarkerlist;	// pointer to face markers
//	int* pmarker = tetgenOut.pointmarkerlist;	// pointer to point markers

    double xyz[6] = { data.xmin, data.xmax, data.ymin, data.ymax, data.zmin, data.zmax };

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        fmarker[i] = 0; // Set default maker value
        // Loop through x, y and z coordinates
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k;	// ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            if (node[l1] == (node[l2] == (node[l3] == xyz[2 * k + 0]))) fmarker[i] = -1 * k; // outer face in negative direction
            if (node[l1] == (node[l2] == (node[l3] == xyz[2 * k + 1]))) fmarker[i] = k;	// outer face in positive direction

        }
        l1 = face[j + 0]; // index of 1st node
        l2 = face[j + 1];	// ..2nd node
        l3 = face[j + 2]; // ..3rd node
//		if ( pmarker[l1] == (pmarker[l2] == (pmarker[l3] == data.type_surf))) fmarker[i] = 100;
    }

    calculate(cmd, &tetgenOut, &tetgenOut);
}

// Function to make the mesh smoother
void Mesher::smooth_mesh(int nr_of_iterations) {
    REAL* node = tetgenOut.pointlist;		// pointer to nodes
    int* face = tetgenOut.trifacelist;		// pointer to triangular faces
    int N = tetgenOut.numberoftrifaces;
    int i, j, k, l1, l2, l3;
    double mean;

    // TODO: exclude boundary nodes from averaging to keep the size of simulation box constant

    // Do several iterations
    for (j = 0; j < nr_of_iterations; ++j) {
        // Loop through the faces
        for (i = 0; i < N; ++i) {
            // Loop through x, y and z coordinates
            for (k = 0; k < 3; ++k) {
                l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
                l2 = 3 * face[j + 1] + k;	// ..2nd node
                l3 = 3 * face[j + 2] + k; // ..3rd node
                mean = (node[l1] + node[l2] + node[l3]) / 3.0;
                node[l1] = node[l2] = node[l3] = mean;
            }
        }
    }
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

// Function to output mesh in .vtk format that can be read by ParaView
void Mesher::output_mesh() {
    calculate("k", &tetgenOut, NULL);
}

// Function to perform tetgen calculation on input and output data
void Mesher::calculate(string cmd, tetgenio* t1, tetgenio* t2) {
    tetgenbeh.parse_commandline(const_cast<char*>(cmd.c_str()));
    tetrahedralize(&tetgenbeh, t1, t2);
}

} /* namespace femocs */
