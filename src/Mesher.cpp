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
#include <vector>
#include <fstream>

#include "Vacuum.h"
#include "Mesh.h"

using namespace std;
namespace femocs {

Mesher::Mesher(string mesher) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
}

// Function to generate simple mesh that consists of one tetrahedron
const shared_ptr<Mesh> Mesher::get_simple_mesh() {
    const int n_nodes = 4;
    const int n_faces = 4;
    const int n_elems = 1;

    shared_ptr<Mesh> new_mesh(new Mesh());

    new_mesh->init_nodes(n_nodes);
    new_mesh->add_node( 1.0, 0.0, 0.7);
    new_mesh->add_node(-1.0, 0.0, 0.7);
    new_mesh->add_node( 0.0, 1.0,-0.7);
    new_mesh->add_node( 0.0,-1.0,-0.7);

    new_mesh->init_faces(n_faces);
    new_mesh->add_face(0,1,3);
    new_mesh->add_face(1,2,3);
    new_mesh->add_face(2,0,3);
    new_mesh->add_face(0,1,2);

    new_mesh->init_elems(n_elems);
    new_mesh->add_elem(0,1,2,3);

    return new_mesh;
}

const shared_ptr<Mesh> Mesher::extract_vacuum_mesh(shared_ptr<Mesh> mesh, const int nmax, const int n_surf, const double latconst, const Femocs::SimuCell* cell, const string cmd) {
    const int n_nodes_per_elem  = 4;

    shared_ptr<Mesh> new_mesh(new Mesh());
    int i, j;
    int n_nodes = mesh->get_n_nodes();
    int n_elems = mesh->get_n_elems();
    int* elems = mesh->get_elems();

    vector<bool> elem_in_vacuum(n_elems);

    int nmax_bulk = nmax - n_surf;
    int n_vacuum, n_not_bottom, n_not_bulk;

    for (i = 0; i < n_elems; ++i) {
        n_vacuum = n_not_bottom = n_not_bulk = 0;

        // Find nodes that are not inside the bulk
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (mesh->get_elem(i,j) >= nmax) n_vacuum++;

        // Find nodes that are not in floating element
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (mesh->get_elem(i,j) >= nmax_bulk) n_not_bulk++;

        // Find elements not belonging to the bottom of simulation cell
        for (j = 0; j < n_nodes_per_elem; ++j) {
            int node = mesh->get_elem(i,j);
            if ( fabs(mesh->get_z(node) - cell->zmin ) > 1.0*latconst ) n_not_bottom++;
        }

        elem_in_vacuum[i] = (n_vacuum > 0)&&(n_not_bulk > 1)&&(n_not_bottom == n_nodes_per_elem);
    }

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    new_mesh->init_nodes(n_nodes);
    new_mesh->copy_nodes(mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    new_mesh->init_elems(n_vacuum_elems);

    for (i = 0; i < n_elems; ++i)
        if(elem_in_vacuum[i]) {
            j = n_nodes_per_elem * i;
            new_mesh->add_elem(elems[j+0], elems[j+1], elems[j+2], elems[j+3]);
        }

    new_mesh->recalc(cmd);
    return new_mesh;
}

// Function to generate union mesh from bulk and vacuum meshes
const shared_ptr<Mesh> Mesher::get_union_mesh(shared_ptr<Mesh> mesh_bulk,
        shared_ptr<Mesh> mesh_volume, const Femocs::SimuCell* cell) {

    shared_ptr<Mesh> new_mesh(new Mesh());
    int i, j, f1, f2, f3;
    int nelems = mesh_volume->get_n_elems();
    int nnodes_volume = mesh_volume->get_n_nodes();
    int nfaces_bulk = mesh_bulk->get_n_facemarkers();
    int nfaces_volume = mesh_volume->get_n_faces();

//    double* nodes = mesh_bulk->getNodes();
    int* faces = mesh_bulk->get_faces();

    vector<bool> is_quality(nfaces_bulk);

    // Mark the faces on the edge of simulation cell
    for (i = 0; i < nfaces_bulk; ++i)
        is_quality[i] = mesh_bulk->get_facemarker(i) == cell->type_surf;

    int nqualityfaces = accumulate(is_quality.begin(), is_quality.end(), 0);

    new_mesh->init_nodes(nnodes_volume);
    new_mesh->init_nodemarkers(nnodes_volume);
    new_mesh->init_faces(nfaces_volume + nqualityfaces);
    new_mesh->init_elems(nelems);

    new_mesh->copy_statistics(mesh_bulk);

    // Copy the nodes of volume mesh without modification
    new_mesh->copy_nodes(mesh_volume);
    new_mesh->copy_nodemarkers(mesh_volume);
    // Copy the faces of volume without modification
    new_mesh->copy_faces(mesh_volume, 0);
    // Copy the elements of volume without modification
    new_mesh->copy_elems(mesh_volume, 0);

    // Copy only the bulk faces that are on the surface of material
    for (i = 0; i < nfaces_bulk; ++i)
        if (is_quality[i]) {
            j = 3 * i;
            f1 = faces[j+0];
            f2 = faces[j+1];
            f3 = faces[j+2];
            new_mesh->add_face(f1, f2, f3);
        }

    return new_mesh;
}

// Function to generate surface faces into available mesh
const shared_ptr<Mesh> Mesher::get_union_mesh_vol2(shared_ptr<Mesh> mesh, const int nmax, const Femocs::SimuCell* cell) {
    const int n_nodes_per_face = 3;
    const int n_nodes_per_elem = 4;

    shared_ptr<Mesh> new_mesh(new Mesh());
    int i, j, n1, n2, n3;
    bool difs[n_nodes_per_elem];
    int n_elems = mesh->get_n_elems();
    int n_nodes = mesh->get_n_nodes();
    int n_faces = mesh->get_n_faces();

    int* elems = mesh->get_elems();
    vector<bool> is_quality(n_elems);

    // Mark the faces on the edge of simulation cell
    for (i = 0; i < n_elems; ++i) {
        for (j = 0; j < n_nodes_per_elem; ++j)
            difs[j] = mesh->get_elem(i,j) < nmax;
        is_quality[i] = difs[0]+difs[1]+difs[2]+difs[3] == 3;
    }

    int n_qualityfaces = accumulate(is_quality.begin(), is_quality.end(), 0);

    new_mesh->init_nodes(n_nodes);
    new_mesh->init_nodemarkers(n_nodes);
    new_mesh->init_elems(n_elems);
    new_mesh->init_faces(n_faces + n_qualityfaces);

    // Copy the nodes of input mesh without modification
    new_mesh->copy_nodes(mesh);
    // Copy the node markers of input mesh without modification
    new_mesh->copy_nodemarkers(mesh);
    // Copy the elements of input mesh without modification
    new_mesh->copy_elems(mesh, 0);
    // Copy the side faces of input mesh
    new_mesh->copy_faces(mesh, 0);

    // Generate the faces that separate material and vacuum
    for (i = 0; i < n_elems; ++i)
        if (is_quality[i]) {
            for (j = 0; j < n_nodes_per_elem; ++j)
                difs[j] = mesh->get_elem(i,j) < nmax;
            j = n_nodes_per_elem * i;

            /* The possible combinations of difs and n1,n2,n3:
             * difs: 1110  1101  1011  0111
             *   n1: elem0 elem0 elem0 elem1
             *   n2: elem1 elem1 elem2 elem2
             *   n3: elem2 elem3 elem3 elem3
             */
            n1 = difs[0] * elems[j+0] + (!difs[0]) * elems[j+1];
            n2 = (difs[0]&difs[1]) * elems[j+1] + (difs[2]&difs[3]) * elems[j+2];
            n3 = (!difs[3]) * elems[j+2] + difs[3] * elems[j+3];
            new_mesh->add_face(n1, n2, n3);
        }

    return new_mesh;
}

// Function to generate mesh from bulk and vacuum atoms
const shared_ptr<Mesh> Mesher::get_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum,
        const string cmd) {
    int i;
    int n_bulk = bulk->get_n_atoms();
    int n_vacuum = vacuum->get_n_atoms();
    shared_ptr<Mesh> new_mesh(new Mesh());

    new_mesh->init_nodes(n_bulk + n_vacuum);
    new_mesh->init_nodemarkers(n_bulk + n_vacuum);
    for (i = 0; i < n_bulk; ++i) {
        new_mesh->add_node(bulk->get_x(i), bulk->get_y(i), bulk->get_z(i));
        new_mesh->add_nodemarker(bulk->get_type(i));
    }
    for (i = 0; i < n_vacuum; ++i) {
        new_mesh->add_node(vacuum->get_x(i), vacuum->get_y(i), vacuum->get_z(i));
        new_mesh->add_nodemarker(vacuum->get_type(i));
    }

    new_mesh->recalc(cmd);
    return new_mesh;
}

// Function to calculate statistics about mesh
void Mesher::calc_statistics(shared_ptr<Mesh> mesh) {
    int N = mesh->get_n_nodes();
    mesh->init_volumes(N);
    mesh->calc_volumes();
    mesh->calc_volume_statistics();
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_faces(shared_ptr<Mesh> mesh, shared_ptr<Surface> surf,
        const Femocs::SimuCell* cell) {
    int i, j, k, l1, l2, l3, m1, m2, facemarker;
    REAL* node = mesh->get_nodes();  // pointer to nodes
    int* face = mesh->get_faces(); // pointer to tetrahedron faces
    int on_xyz_surface;

    const int ncoords = 3; // nr of coordinates
    const int nnodes = 3;  // nr of nodes per face
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);

    const double xyz[6]  = { cell->xmin, cell->xmax, cell->ymin, cell->ymax, cell->zmin,
            cell->zmaxbox };
    const int markers[6] = { cell->type_xmin, cell->type_xmax, cell->type_ymin, cell->type_ymax,
            cell->type_zmin, cell->type_zmax };

    // Loop through the faces
    for (i = 0; i < N; ++i) {
        j = nnodes * i;
        facemarker = cell->type_none; // Set default maker value
        on_xyz_surface = 0;
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
            else if (on_face(node[l1], node[l2], node[l3], xyz[m2]))
                facemarker = markers[m2];

            else if (on_surface(node[l1], node[l2], node[l3], surf, k))
                on_xyz_surface++;

        }
        if(on_xyz_surface == ncoords) facemarker = cell->type_surf;
        mesh->add_facemarker(facemarker);
    }
}

// Function to mark the outer edges of simulation cell
void Mesher::mark_elems(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell) {
    int i, j, k, l1, l2, l3, l4, m1, m2, emarker;
    int N = mesh->get_n_elems();
    const int ncoords = 3; // nr of coordinates
    const int nnodes = 4;  // nr of nodes per element

    REAL* node = mesh->get_nodes();       // pointer to nodes
    int* elem = mesh->get_elems();  // pointer to tetrahedrons
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

// Function to mark the outer edges of simulation cell
void Mesher::mark_elems_byvol(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell) {
    int N = mesh->get_n_elems();

    mesh->init_volumes(N);
    mesh->calc_volumes();
    mesh->init_elemmarkers(N);

    const double Vmin = mesh->stat.Vmin;
    const double Vmax = mesh->stat.Vmax;

    for (int i = 0; i < N; ++i) {
        if(mesh->get_volume(i) >= Vmin && mesh->get_volume(i) <= Vmax)
            mesh->add_elemmarker(cell->type_bulk);
        else
            mesh->add_elemmarker(cell->type_vacuum);
    }
}

// Mark elements by comparing the centres of elements with the ones in bulk mesh
void Mesher::mark_elems_bycentre(shared_ptr<Mesh> mesh, shared_ptr<Mesh> bulk_mesh, const Femocs::SimuCell* cell) {
    int N = mesh->get_n_elems();

    mesh->calc_centres();
    bulk_mesh->calc_centres();

    mesh->init_elemmarkers(N);

    for (int i = 0; i < N; ++i)
        if( mesh->centre_found_in(i, bulk_mesh) )
            mesh->add_elemmarker(cell->type_bulk);
        else
            mesh->add_elemmarker(cell->type_vacuum);

}

/*
    The idea is, that when tetgen adds nodes to the mesh, it add them to the end of nodelist
    In that ways the first nodes will always be the bulk (or surface) nodes and there is
    no need to search them later on. Those nodes can be used to find the faces on the surface.
*/
void Mesher::mark_elems_bynode(shared_ptr<Mesh> mesh, const int nmax, const Femocs::SimuCell* cell) {
    int N = mesh->get_n_elems();
    mesh->init_elemmarkers(N);
    bool difs[4];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 4; ++j)
            difs[j] = mesh->get_elem(i,j) < nmax;
    
        if( difs[0]+difs[1]+difs[2]+difs[3] == 3 ) {
            //mesh->add_elemmarker(cell->type_surf);
            if     (difs[0] & difs[3])
                mesh->add_elemmarker(1);
            else if(difs[0] & (!difs[3]))
                mesh->add_elemmarker(2);
            else if((!difs[0]) & difs[3])
                mesh->add_elemmarker(2);
            else if((!difs[0]) & (!difs[3]))
                mesh->add_elemmarker(4);
        }
        else
            //mesh->add_elemmarker(cell->type_none);
            mesh->add_elemmarker(0);
    }
}

// Mark elements by comparing the centres of elements with the ones in bulk mesh
void Mesher::mark_faces_bynode(shared_ptr<Mesh> mesh, const int nmax, const Femocs::SimuCell* cell) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);
    const int n_nodes_per_face = 3;
    bool difs[n_nodes_per_face];
    
    int* face = mesh->get_faces();

    vector<int> face_indxs(N);

    for(int i = 0; i < n_nodes_per_face*N; ++i)
        face_indxs[face[i]] = i;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < n_nodes_per_face; ++j)
            difs[j] = mesh->get_face(i,j) < nmax;

        if( difs[0]+difs[1]+difs[2] == 3 )
            //mesh->add_elemmarker(cell->type_bulk);
            mesh->add_facemarker(cell->type_surf);
        else
            //mesh->add_elemmarker(cell->type_vacuum);
            mesh->add_facemarker(cell->type_none);
    }
}

// Mark elements by comparing the centres of elements with the ones in bulk mesh
void Mesher::mark_faces_bysequence(shared_ptr<Mesh> mesh, const int nmax, const Femocs::SimuCell* cell) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);
    const int n_nodes_per_face = 3;
    bool difs[n_nodes_per_face];

    int* face = mesh->get_faces();

    vector<int> face_indxs(N);

    for(int i = 0; i < n_nodes_per_face*N; ++i)
        face_indxs[face[i]] = (int) i/n_nodes_per_face;

    for (int i = 0; i < N; ++i)
        mesh->add_facemarker(cell->type_none);

    for (int i = (N-1); i >= 0; --i) {
        for (int j = 0; j < n_nodes_per_face; ++j)
            difs[j] = mesh->get_face(i,j) < nmax;
    
        if( difs[0]+difs[1]+difs[2] == 3 )
            mesh->set_facemarker(face_indxs[i], cell->type_surf);
    }
}

// Determine whether the nodes of two faces are overlapping
bool Mesher::on_face(const double f1n1, const double f1n2, const double f1n3, const double f2) {
    double eps = 0.2;
    bool dif1 = (fabs(f1n1 - f2) < eps);
    bool dif2 = (fabs(f1n2 - f2) < eps);
    bool dif3 = (fabs(f1n3 - f2) < eps);
    return dif1 && dif2 && dif3;
}

inline bool point_on_surface(vector<double>* nodes, double node, const int N) {
    const double eps = 0.5;
    for (int i = 0; i < N; ++i)
        if ( fabs((*nodes)[i] - node) < eps )
            return true;
    return false;
}

// Determine whether the bulk mesh face is on the surface
bool Mesher::on_surface(const double f1n1, const double f1n2, const double f1n3,
        shared_ptr<Surface> surf, const int xn) {

    vector<double>* nodes;
    int N = surf->get_n_atoms();

    if (xn == 0)
        nodes = surf->get_xs();
    else if (xn == 1)
        nodes = surf->get_ys();
    else if (xn == 2)
        nodes = surf->get_zs();
    else
        return false;

    bool dif1 = point_on_surface(nodes, f1n1, N);
    bool dif2 = point_on_surface(nodes, f1n2, N);
    bool dif3 = point_on_surface(nodes, f1n3, N);

    return dif1 & dif2 & dif3;
}

// Function to remove too big tetrahedra from the mesh
void Mesher::clean_elems(shared_ptr<Mesh> mesh, const double rmax, const string cmd) {
    REAL* node = mesh->get_nodes();      // pointer to nodes
    int* elem = mesh->get_elems(); // pointer to tetrahedron corners

    double dx, r12, r13, r14, r23, r24, r34;
    int i, j, k, l1, l2, l3, l4;

    int N = mesh->get_n_elems();
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
    temp_mesh->init_elems(mesh->get_n_elems());
    temp_mesh->copy_elems(mesh, 0);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_elems(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert tetrahedra with suitable quality
    update_list(mesh->get_elems(), temp_mesh->get_elems(), isQuality, 4);

    mesh->recalc(cmd);
}

// Function to remove too big faces from the mesh
void Mesher::clean_faces(shared_ptr<Mesh> mesh, const double rmax, const string cmd) {
    REAL* node = mesh->get_nodes();		// pointer to nodes
    int* face = mesh->get_faces();		// pointer to triangular faces

    double dx, r12, r13, r23;
    int i, j, k, l1, l2, l3;
    int N = mesh->get_n_faces();
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
    temp_mesh->init_faces(mesh->get_n_faces());
    temp_mesh->copy_faces(mesh, 0);

    // Initialise and fill the input mesh with cleaned elements
    mesh->init_faces(accumulate(isQuality.begin(), isQuality.end(), 0));

    // Insert faces with suitable quality
    update_list(mesh->get_faces(), temp_mesh->get_faces(), isQuality, 3);

    mesh->recalc(cmd);
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
