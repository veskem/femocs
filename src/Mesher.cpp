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
#include "Mesh.h"

using namespace std;
namespace femocs {

Mesher::Mesher(string mesher) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
}

// Function to generate union mesh from bulk and vacuum meshes
const shared_ptr<Mesh> Mesher::get_union_mesh(shared_ptr<Mesh> mesh_bulk,
        shared_ptr<Mesh> mesh_volume, const Femocs::SimuCell* cell, const string cmd) {

    shared_ptr<Mesh> new_mesh(new Mesh());
    int i, j, f1, f2, f3;

    int nelems = mesh_volume->getNelems();
    int nnodes = mesh_volume->getNnodes();
    int nfaces_bulk = mesh_bulk->getNfacemarkers();
    int nfaces_volume = mesh_volume->getNfaces();

    vector<bool> is_quality(nfaces_bulk);

    // Mark the faces on the edge of simulation cell
    for (i = 0; i < nfaces_bulk; ++i)
        is_quality[i] = mesh_bulk->getFacemarker(i) == cell->type_surf;

    int nfaces = nfaces_volume + accumulate(is_quality.begin(), is_quality.end(), 0);

    // Copy the points of volume mesh without modification
    new_mesh->init_nodes(nnodes);
    new_mesh->copy_nodes(mesh_volume);

    // Copy the elements of bulk without modification
    new_mesh->init_elems(nelems);
    new_mesh->copy_elems(mesh_volume);

    // Copy the faces of bulk without modification
    new_mesh->init_faces(nfaces);
    new_mesh->copy_faces(mesh_volume);

    // Copy only the bulk faces that are not on the edge of simulation cell
    for (i = 0; i < nfaces_bulk; ++i)
        if (is_quality[i]) {
            j = 3 * i;
            f1 = mesh_bulk->getFaces()[j + 0];
            f2 = mesh_bulk->getFaces()[j + 1];
            f3 = mesh_bulk->getFaces()[j + 2];
            new_mesh->add_face(f1, f2, f3);
        }
    new_mesh->calculate(cmd);
    return new_mesh;
}

// Function to generate mesh from bulk and vacuum atoms
const shared_ptr<Mesh> Mesher::get_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum,
        const string cmd) {
    int i;
    int nbulk = bulk->getN();
    int nvacuum = vacuum->getN();
    shared_ptr<Mesh> new_mesh(new Mesh());

    new_mesh->init_nodes(nbulk + nvacuum);
    for (i = 0; i < nbulk; ++i)
        new_mesh->add_node(bulk->getX(i), bulk->getY(i), bulk->getZ(i));
    for (i = 0; i < nvacuum; ++i)
        new_mesh->add_node(vacuum->getX(i), vacuum->getY(i), vacuum->getZ(i));

    new_mesh->calculate(cmd);
    return new_mesh;
}

// Function to generate mesh from bulk atoms
const shared_ptr<Mesh> Mesher::get_bulk_mesh(shared_ptr<Surface> bulk, const string cmd) {
    int i;
    int nbulk = bulk->getN();
    shared_ptr<Mesh> new_mesh(new Mesh());

    new_mesh->init_nodes(nbulk);
    for (i = 0; i < nbulk; ++i)
        new_mesh->add_node(bulk->getX(i), bulk->getY(i), bulk->getZ(i));

    new_mesh->calculate(cmd);
    return new_mesh;
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_faces(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell) {
    int i, j, k, l1, l2, l3, m1, m2, facemarker;
    REAL* node = mesh->getNodes();  // pointer to nodes
    int* face = mesh->getFaces(); // pointer to tetrahedron faces

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 3;  // nr of nodes per face
    int N = mesh->getNfaces();
    mesh->init_facemarkers(N);

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };
    const int markers[6] = { cell->type_xmin, cell->type_xmax, cell->type_ymin, cell->type_ymax,
            cell->type_zmin, cell->type_zmax };

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        facemarker = cell->type_surf; // Set default maker value

        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * face[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * face[j + 1] + k; // ..2nd node
            l3 = ncoords * face[j + 2] + k; // ..3rd node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // outer face in negative direction
            if (on_face(node[l1], node[l2], node[l3], xyz[m1]))
                facemarker = markers[m1];
            // outer face in positive direction
            else if (on_face(node[l1], node[l2], node[l3], xyz[m2])) facemarker = markers[m2];
        }
        mesh->add_facemarker(facemarker);
    }
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_elems(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell) {
    int i, j, k, l1, l2, l3, l4, m1, m2, emarker;
    int N = mesh->getNelems();
    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element

    REAL* node = mesh->getNodes();       // pointer to nodes
    int* elem = mesh->getElems();  // pointer to tetrahedrons
    mesh->init_elemmarkers(N);

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };
    const int markers[6] = { cell->type_xmin, cell->type_xmax, cell->type_ymin, cell->type_ymax,
            cell->type_zmin, cell->type_zmax };

    // Loop through the elements
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        emarker = cell->type_none; // Set default maker value
        // Loop through x, y and z coordinates
        for (k = 0; k < ncoords; ++k) {
            l1 = ncoords * elem[j + 0] + k; // index of x, y or z coordinate of 1st node
            l2 = ncoords * elem[j + 1] + k; // ..2nd node
            l3 = ncoords * elem[j + 2] + k; // ..3rd node
            l4 = ncoords * elem[j + 3] + k; // ..4th node

            m1 = 2 * k + 0; // index for min x, y or z-coordinate
            m2 = 2 * k + 1; //           max x, y or z-coordinate

            // element on negative direction edge
            if (on_face(node[l1], node[l2], node[l3], xyz[m1])
                    || on_face(node[l1], node[l2], node[l4], xyz[m1])
                    || on_face(node[l1], node[l3], node[l4], xyz[m1])
                    || on_face(node[l2], node[l3], node[l4], xyz[m1])) emarker = markers[m1];

            // element on positive direction edge
            if (on_face(node[l1], node[l2], node[l3], xyz[m2])
                    || on_face(node[l1], node[l2], node[l4], xyz[m2])
                    || on_face(node[l1], node[l3], node[l4], xyz[m2])
                    || on_face(node[l2], node[l3], node[l4], xyz[m2])) emarker = markers[m2];
        }
        mesh->add_elemmarker(emarker);
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

// Determine whether the bulk mesh face is on the surface
bool Mesher::on_surface(const double f1n1, const double f1n2, const double f1n3,
        const double zmin) {
    return (f1n1 >= zmin) && (f1n2 >= zmin) && (f1n3 >= zmin);
}

// Function to remove too big tetrahedra from the mesh
void Mesher::clean_elems(shared_ptr<Mesh> mesh, const double rmax, const string cmd) {
    REAL* node = mesh->getNodes();      // pointer to nodes
    int* elem = mesh->getElems(); // pointer to tetrahedron corners

    double dx, r12, r13, r14, r23, r24, r34;
    int i, j, k, l1, l2, l3, l4;

    int N = mesh->getNelems();
    vector<bool> isQuality(N);

    // Loop through the tetrahedra
    for (i = 0; i < N; ++i) {
        j = 4 * i;
        // Loop through x, y and z coordinates
        r12 = r13 = r14 = r23 = r24 = r34 = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * elem[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * elem[j + 1] + k;   // ..2nd node
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

    // Make temporary copy from the input mesh
    shared_ptr<Mesh> temp_mesh(new Mesh());
    temp_mesh->init_elems(mesh->getNelems());
    temp_mesh->copy_elems(mesh);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_elems(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert tetrahedra with suitable quality
    update_list(mesh->getElems(), temp_mesh->getElems(), isQuality, 4);

    mesh->calculate(cmd);
}

// Function to remove too big faces from the mesh
void Mesher::clean_faces(shared_ptr<Mesh> mesh, const double rmax, const string cmd) {
    REAL* node = mesh->getNodes();		// pointer to nodes
    int* face = mesh->getFaces();		// pointer to triangular faces

    double dx, r12, r13, r23;
    int i, j, k, l1, l2, l3;
    int N = mesh->getNfaces();
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

    // Make temporary copy from the input mesh
    shared_ptr<Mesh> temp_mesh(new Mesh());
    temp_mesh->init_faces(mesh->getNfaces());
    temp_mesh->copy_faces(mesh);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_faces(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert faces with suitable quality
    update_list(mesh->getFaces(), temp_mesh->getFaces(), isQuality, 3);

    mesh->calculate(cmd);
}

// Function to remove too big faces from the mesh
// SHOULD BUT FOR SOME REASON DOESN'T REMOVE TO FACES ON SIMU CELL EDGES
void Mesher::clean_facets(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell, const string cmd) {
    REAL* node = mesh->getNodes();      // pointer to nodes
    int* face = mesh->getFaces();     // pointer to triangular faces

    int i, j, k, l1, l2, l3;
    bool on_simucell_edge;
    int N = mesh->getNfaces();
    vector<bool> isQuality(N);

    const double xyz[6] = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = 3 * i;
        // Loop through x, y and z coordinates
        on_simucell_edge = 0;
        for (k = 0; k < 3; ++k) {
            l1 = 3 * face[j + 0] + k; // index of x,y or z coordinate of 1st node
            l2 = 3 * face[j + 1] + k; // ..2nd node
            l3 = 3 * face[j + 2] + k; // ..3rd node

            // Detect whether the face is on the edge of simulation box
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2 * k + 0]);
            on_simucell_edge |= on_face(node[l1], node[l2], node[l3], xyz[2 * k + 1]);
        }
        // Keep only the faces with appropriate quality
        isQuality[i] = (!on_simucell_edge);
    }

    // Make temporary copy from the input mesh
    shared_ptr<Mesh> temp_mesh(new Mesh());
    temp_mesh->init_faces(mesh->getNfaces());
    temp_mesh->copy_faces(mesh);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_faces(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert faces with suitable quality
    update_list(mesh->getFaces(), temp_mesh->getFaces(), isQuality, 3);
    mesh->calculate(cmd);
}

// Function to remove the objects from element, face or edge list
void Mesher::update_list(int* new_list, int* old_list, vector<bool> is_quality, int M) {
    int i, j, k;
    int N = is_quality.size();
    j = 0;
    // loop through the old array and move quality objects to the left
    for (i = 0; i < N; ++i) {
        if (is_quality[i]) {
            for (k = 0; k < M; ++k)
                new_list[j + k] = old_list[M * i + k];
            j += M;
        }
    }
}

} /* namespace femocs */
