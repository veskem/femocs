/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>

using namespace std;
namespace femocs {

// Function to calculate the indices of sorted array
template<class Vals>
void get_sort_permutation(const Vals& values, std::vector<int>& v) {
    int size = values.size();
    v.clear();
    v.reserve(size);
    for (int i = 0; i < size; ++i)
        v.push_back(i);

    sort(v.begin(), v.end(), [&values](int a, int b) -> bool {
        return values[a] < values[b];
    });
}

Mesher::Mesher(string mesher) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
}

// Function to generate mesh between surface atoms
void Mesher::generate_mesh(shared_ptr<Surface> bulk, shared_ptr<Surface> surf, string cmd) {
//    generate_surface_mesh(surf, &tetgenIn);
    generate_volume_mesh(bulk, surf, &tetgenIn);
    calculate(cmd, &tetgenIn, &tetgenOut);
}

// Function to sort the surface atoms according to their radial distance on xy-plane
void Mesher::sort_surface(vector<int>& permutation_indxs, shared_ptr<Surface> surf) {
    int i;
    int N = surf->getN();
    double x, y;
    vector<double> r2;
    r2.reserve(N);
    for (i = 0; i < N; ++i) {
        x = surf->getX(i);
        y = surf->getY(i);
        r2.push_back(x*x + y*y);
    }
    get_sort_permutation(r2, permutation_indxs);

//    for (i = 0; i < N; ++i)
//        cout << permutation_indxs[i] << ": " << r2[permutation_indxs[i]] << endl;
}
/*
 * TODO: add bulk material atoms to the mesh and avoid using surface atoms twice when
 * generating vacuum virtual atoms as surface atoms elevated in z direction
 * Also replace the elevation factor from 2^i to double^i

 * To remove every second atom from surface, they must be sorted first.
 * It can be done by first calculating the atomic radial coordinates in
 * xy-plane, then finding sort indices for them (like described here)
 * http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 * and finally sorting(or just using sort indices for) the x,y,z,type... vectors.
 */

// Function to generate volumetric mesh between surface atoms and vacuum
void Mesher::generate_volume_mesh(shared_ptr<Surface> bulk, shared_ptr<Surface> surf, tetgenio* tetIO) {
    int i, j, k1, k2;
    int nlayers = 1;            // number of node layers
    int N_bulk = bulk->getN();  // number of atoms in bulk material
    int N_surf = surf->getN();  // number of atoms in surface

    // sum of Nsurf/2^1 + Nsurf/2^2 +...+ Nsurf/2^nlayers = Nsurf*( 1 - 1/2^nlayers )/(1 - 1/2)
    int M = ceil(4*nlayers + N_bulk + N_surf*(1.0 - 1.0/(1 << nlayers))); // total number of nodes

    tetIO->numberofpoints = M;
    tetIO->pointlist = new REAL[3 * M];
    tetIO->pointmarkerlist = new int[M];

    // TODO!!! GET THOSE CORNERS AUTOMATICALLY!!!
    const double pos_corn = 74.99;
    const double neg_corn = -74.99;
    const double zmax = 100;
    const double dz = 15;
    double znew = zmax;

    vector<int> iperm;
    sort_surface(iperm, surf);
    k1 = 0; k2 = 0;

    // Add bulk atoms
    k1 = 0;
    for (i = 0; i < N_bulk; ++i) {
        add_point(bulk->getX(i), bulk->getY(i), bulk->getZ(i), bulk->getType(i), k1, tetIO);
        k1 += 3;
    }

    // Add every 2nd atom to the next layer as compared to previous one
    for (i = 1; i <= nlayers; ++i) {
        for (j = 0; j < N_surf / (1 << i); ++j) {
            k2 = iperm[j * (1 << i)];
            //k2 = j * (1 << i);

            //cout << k1 << " " << k2 << " " << znew << endl;
            if (k1 >= 3 * M) {
                cout << "k1 out of limits!, i,j=" << i << "," << j << endl; exit(EXIT_FAILURE);
            }
            if (k2 >= N_surf) {
                cout << "k2 out of limits!" << endl; exit(EXIT_FAILURE);
            }

            znew = dz * i + surf->getZ(k2);
//            znew = surf->getZ(k2) + pow(dz, i);
            if (znew > zmax) znew = zmax;

            add_point(surf->getX(k2), surf->getY(k2), znew, surf->getType(k2), k1, tetIO);
            k1 += 3;
        }

        add_point(pos_corn, pos_corn, znew, 1, k1, tetIO); k1 += 3;
        add_point(pos_corn, neg_corn, znew, 1, k1, tetIO); k1 += 3;
        add_point(neg_corn, pos_corn, znew, 1, k1, tetIO); k1 += 3;
        add_point(neg_corn, neg_corn, znew, 1, k1, tetIO); k1 += 3;
    }

    // Compensate the mismatch between allocated array size and actually inserted points
    // Mismatch is caused by the rounding error in double->int conversion
    while(k1 / 3 < M) {
        add_point(surf->getX(0), surf->getY(0), surf->getZ(0), surf->getType(0), k1, tetIO);
        k1 += 3;
    }
}

// Function to add point with coordinates and maker to tetgen list
void Mesher::add_point(const double x, const double y, const double z, const int pmarker, const int i, tetgenio* tetIO) {
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

// Function to mark the outer edges of simulation cell
void Mesher::mark_mesh(AtomReader::Data data, string cmd) {
    int i, j, k, l1, l2, l3, m1, m2;
    int N = tetgenOut.numberoftrifaces;

    REAL* node = tetgenOut.pointlist;	// pointer to nodes
    int* face = tetgenOut.trifacelist;	// pointer to tetrahedron faces

    tetgenOut.trifacemarkerlist = new int[N];
    int* fmarker = tetgenOut.trifacemarkerlist; // pointer to face markers
    tetgenIn.pointmarkerlist = new int[tetgenOut.numberofpoints];
//    int* pmarker = tetgenOut.pointmarkerlist;   // pointer to point markers

    double eps = 0.1;

    double xyz[6] = { data.xmin, data.xmax, data.ymin, data.ymax, data.zmin, data.zmax };
    bool dif[6];

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        fmarker[i] = 66; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k;	// ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            m1 = 2 * k + 0;
            m2 = 2 * k + 1;

            dif[0] = (fabs(node[l1] - xyz[m1]) < eps);
            dif[1] = (fabs(node[l2] - xyz[m1]) < eps);
            dif[2] = (fabs(node[l3] - xyz[m1]) < eps);

            dif[3] = (fabs(node[l1] - xyz[m2]) < eps);
            dif[4] = (fabs(node[l2] - xyz[m2]) < eps);
            dif[5] = (fabs(node[l3] - xyz[m2]) < eps);

            if (dif[0] && dif[1] && dif[2]) {
                //      cout << node[l1] << " " << node[l2] << " " << node[l3] << endl;
                fmarker[i] = 10 * (k + 1); // outer face in negative direction
//                pmarker[face[j+0]] = 10*(k+1);
//                pmarker[face[j+1]] = 10*(k+1);
//                pmarker[face[j+2]] = 10*(k+1);
            }

            if (dif[3] && dif[4] && dif[5]) {
                fmarker[i] = k + 1;   // outer face in positive direction
//                pmarker[face[j+0]] = (k+1);
//                pmarker[face[j+1]] = (k+1);
//                pmarker[face[j+2]] = (k+1);
            }

        }
//        l1 = face[j + 0]; // index of 1st node
//        l2 = face[j + 1];	// ..2nd node
//        l3 = face[j + 2]; // ..3rd node
//        if ( pmarker[l1] == (pmarker[l2] == (pmarker[l3] == data.type_surf))) fmarker[i] = 100;
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
