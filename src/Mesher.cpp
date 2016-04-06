/*
 * Mesher.cpp
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#include "Mesher.h"

using namespace std;
namespace femocs {

Mesher::Mesher(string mesher, const double latconst) {
    if (mesher != "tetgen") cout << "Unknown mesher: " + mesher << endl;
    this->latconst = latconst;
}

// Function to generate simple mesh that consists of one tetrahedron
const void Mesher::get_test_mesh(Mesh* new_mesh) {
    const int n_nodes = 4;
    const int n_faces = 4;
    const int n_elems = 1;

    new_mesh->init_nodes(n_nodes);
    new_mesh->add_node(1.0, 0.0, 0.7);
    new_mesh->add_node(-1.0, 0.0, 0.7);
    new_mesh->add_node(0.0, 1.0, -0.7);
    new_mesh->add_node(0.0, -1.0, -0.7);

    new_mesh->init_faces(n_faces);
    new_mesh->add_face(0, 1, 3);
    new_mesh->add_face(1, 2, 3);
    new_mesh->add_face(2, 0, 3);
    new_mesh->add_face(0, 1, 2);

    new_mesh->init_elems(n_elems);
    new_mesh->add_elem(0, 1, 2, 3);
}

// Function to generate mesh from bulk and vacuum atoms
const void Mesher::get_volume_mesh(Mesh* new_mesh, Bulk* bulk_mesh, Vacuum* vacuum_mesh,
        const string cmd) {
    int i;
    int n_bulk = bulk_mesh->get_n_atoms();
    int n_vacuum = vacuum_mesh->get_n_atoms();

    new_mesh->init_nodes(n_bulk + n_vacuum);

    for (i = 0; i < n_bulk; ++i)
        new_mesh->add_node(bulk_mesh->get_x(i), bulk_mesh->get_y(i), bulk_mesh->get_z(i));

    for (i = 0; i < n_vacuum; ++i)
        new_mesh->add_node(vacuum_mesh->get_x(i), vacuum_mesh->get_y(i), vacuum_mesh->get_z(i));

    new_mesh->recalc(cmd);
}

const void Mesher::extract_mesh(vector<bool>* is_vacuum, Mesh* big_mesh, const int nmax,
        const int n_surf, const double zmin) {
    int i, j;
    int n_elems = big_mesh->get_n_elems();

    is_vacuum->resize(n_elems);

    int nmax_bulk = nmax - n_surf;
    int n_vacuum, n_not_bottom, n_not_bulk;

    for (i = 0; i < n_elems; ++i) {
        n_vacuum = n_not_bottom = n_not_bulk = 0;

        // Find nodes that are not inside the bulk
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (big_mesh->get_elem(i, j) >= nmax) n_vacuum++;

        // Find nodes that are not in floating element
        for (j = 0; j < n_nodes_per_elem; ++j)
            if (big_mesh->get_elem(i, j) >= nmax_bulk) n_not_bulk++;

        // Find elements not belonging to the bottom of simulation cell
        for (j = 0; j < n_nodes_per_elem; ++j) {
            int node = big_mesh->get_elem(i, j);
            if (fabs(big_mesh->get_z(node) - zmin) > 1.0 * latconst) n_not_bottom++;
        }

        (*is_vacuum)[i] = (n_vacuum > 0) && (n_not_bulk > 1) && (n_not_bottom == n_nodes_per_elem);
    }
}

const void Mesher::separate_vacuum_mesh(Mesh* vacuum_mesh, Mesh* big_mesh, const int nmax,
        const int n_surf, const double zmin, const string cmd) {
    int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, nmax, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    vacuum_mesh->init_nodes(big_mesh->get_n_nodes());
    vacuum_mesh->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    vacuum_mesh->init_elems(n_vacuum_elems);

    for (int i = 0; i < big_mesh->get_n_elems(); ++i)
        if (elem_in_vacuum[i]) {
            int j = n_nodes_per_elem * i;
            vacuum_mesh->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    vacuum_mesh->recalc(cmd);
}

const void Mesher::separate_bulk_mesh(Mesh* bulk, Mesh* big_mesh, const int nmax, const int n_surf,
        const double zmin, const string cmd) {
    int* elems = big_mesh->get_elems();
    int n_elems = big_mesh->get_n_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, nmax, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    bulk->init_nodes(big_mesh->get_n_nodes());
    bulk->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    bulk->init_elems(n_elems - n_vacuum_elems);

    for (int i = 0; i < n_elems; ++i)
        if (!elem_in_vacuum[i]) {
            int j = n_nodes_per_elem * i;
            bulk->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    bulk->recalc(cmd);
}

const void Mesher::separate_meshes(Mesh* bulk, Mesh* vacuum, Mesh* big_mesh, const int nmax,
        const int n_surf, const double zmin, const string cmd) {

    int i, j;
    int n_nodes = big_mesh->get_n_nodes();
    int n_elems = big_mesh->get_n_elems();
    int* elems = big_mesh->get_elems();

    vector<bool> elem_in_vacuum;
    extract_mesh(&elem_in_vacuum, big_mesh, nmax, n_surf, zmin);

    int n_vacuum_elems = accumulate(elem_in_vacuum.begin(), elem_in_vacuum.end(), 0);

    // Copy the nodes from input mesh without modification
    vacuum->init_nodes(n_nodes);
    vacuum->copy_nodes(big_mesh);

    bulk->init_nodes(n_nodes);
    bulk->copy_nodes(big_mesh);

    // Copy only the elements from input mesh that have the node outside the bulk nodes
    vacuum->init_elems(n_vacuum_elems);
    bulk->init_elems(n_elems - n_vacuum_elems);

    for (i = 0; i < n_elems; ++i)
        if (elem_in_vacuum[i]) {
            j = n_nodes_per_elem * i;
            vacuum->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        } else {
            j = n_nodes_per_elem * i;
            bulk->add_elem(elems[j + 0], elems[j + 1], elems[j + 2], elems[j + 3]);
        }

    vacuum->recalc(cmd);
    bulk->recalc(cmd);
}

// Function to generate surface faces into available mesh
const void Mesher::generate_surf_faces(Mesh* mesh, const int nmax) {
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
            difs[j] = mesh->get_elem(i, j) < nmax;
        is_quality[i] = difs[0] + difs[1] + difs[2] + difs[3] == 3;
    }

    int n_qualityfaces = accumulate(is_quality.begin(), is_quality.end(), 0);

    // Make temporary copy of available faces
    Mesh temp_mesh;
    temp_mesh.init_faces(n_faces);
    temp_mesh.copy_faces(mesh, 0);

    // Make face list size bigger
    mesh->init_faces(n_faces + n_qualityfaces);
    // Copy back already available faces
    mesh->copy_faces(&temp_mesh, 0);

    // Generate the faces that separate material and vacuum
    for (i = 0; i < n_elems; ++i)
        if (is_quality[i]) {
            for (j = 0; j < n_nodes_per_elem; ++j)
                difs[j] = mesh->get_elem(i, j) < nmax;
            j = n_nodes_per_elem * i;

            /* The possible combinations of difs and n1,n2,n3:
             * difs: 1110  1101  1011  0111
             *   n1: elem0 elem0 elem0 elem1
             *   n2: elem1 elem1 elem2 elem2
             *   n3: elem2 elem3 elem3 elem3
             */
            n1 = difs[0] * elems[j + 0] + (!difs[0]) * elems[j + 1];
            n2 = (difs[0] & difs[1]) * elems[j + 1] + (difs[2] & difs[3]) * elems[j + 2];
            n3 = (!difs[3]) * elems[j + 2] + difs[3] * elems[j + 3];
            mesh->add_face(n1, n2, n3);
        }
}

/*
 The idea is, that when tetgen adds nodes to the mesh, it adds them to the end of node list
 In that ways the first nodes will always be the bulk (or surface) nodes and there is
 no need to search them later on. Those nodes can be used to find the faces on the surface.
 */
void Mesher::mark_elems_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types) {
    int N = mesh->get_n_elems();
    mesh->init_elemmarkers(N);
    bool difs[4];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 4; ++j)
            difs[j] = mesh->get_elem(i, j) < nmax;

        if (difs[0] + difs[1] + difs[2] + difs[3] == 3) {
            if (difs[0] & difs[3])
                mesh->add_elemmarker(1);
            else if (difs[0] & (!difs[3]))
                mesh->add_elemmarker(2);
            else if ((!difs[0]) & difs[3])
                mesh->add_elemmarker(2);
            else if ((!difs[0]) & (!difs[3])) mesh->add_elemmarker(4);
        } else
            mesh->add_elemmarker(types->type_none);
    }
}

// Mark elements by comparing the centres of elements with the ones in bulk mesh
void Mesher::mark_faces_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);
    bool difs[n_nodes_per_face];

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < n_nodes_per_face; ++j)
            difs[j] = mesh->get_face(i, j) < nmax;

        if (difs[0] + difs[1] + difs[2] == 3)
            mesh->add_facemarker(types->type_surf);
        else
            mesh->add_facemarker(types->type_none);
    }
}

// Function to determine whether the center of face is on the boundary of simulation cell or not
inline bool on_boundary(const double face, const double face_max) {
    const double eps = 0.1;
    return fabs(face - face_max) < eps;
}

// Mark the boundary faces of mesh
void Mesher::mark_faces(Mesh* mesh, const AtomReader::Sizes* sizes, const AtomReader::Types* types) {
    int N = mesh->get_n_faces();
    mesh->init_facemarkers(N);

    for (int i = 0; i < N; ++i) {
        if (on_boundary(mesh->get_face_centre(i, 0), sizes->xmin))
            mesh->add_facemarker(types->type_xmin);
        else if (on_boundary(mesh->get_face_centre(i, 0), sizes->xmax))
            mesh->add_facemarker(types->type_xmax);
        else if (on_boundary(mesh->get_face_centre(i, 1), sizes->ymin))
            mesh->add_facemarker(types->type_ymin);
        else if (on_boundary(mesh->get_face_centre(i, 1), sizes->ymax))
            mesh->add_facemarker(types->type_ymax);
        else if (on_boundary(mesh->get_face_centre(i, 2), sizes->zminbox))
            mesh->add_facemarker(types->type_zmin);
        else if (on_boundary(mesh->get_face_centre(i, 2), sizes->zmaxbox))
            mesh->add_facemarker(types->type_zmax);
        else
            mesh->add_facemarker(types->type_surf);
    }
}

// Function to remove too big tetrahedra from the mesh
void Mesher::clean_elems(Mesh* mesh, const double rmax, const string cmd) {
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
void Mesher::clean_faces(Mesh* mesh, const double rmax, const string cmd) {
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
