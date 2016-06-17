/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include <fstream>

using namespace std;
namespace femocs {

// SolutionReader constructor
SolutionReader::SolutionReader() :
        longrange_efield(0) {
    reserve(0);
}

const void SolutionReader::extract_solution(DealII* fem, Medium &medium) {
    dealii::Tensor<1, DIM> ef;
    double pot;

    this->fem = fem;
    longrange_efield = fem->get_efield();
    int n_nodes = medium.get_n_atoms();
    reserve(n_nodes);

    vector<int> medium2node_map = get_medium2node_map(medium);
    vector<int> node2elem_map = get_node2elem_map();
    vector<int> node2vert_map = get_node2vert_map();

    for (int node = 0; node < n_nodes; ++node) {
        int n = medium2node_map[node];
        if (n >= 0) {
            ef = fem->get_elfield_at_node(node2elem_map[n], node2vert_map[n]);
            pot = fem->get_potential_at_node(node2elem_map[n], node2vert_map[n]);
            add_atom(medium.get_id(node), medium.get_point(node), medium.get_coordination(node));
        } else {
            ef[0] = ef[1] = ef[2] = pot = 1e20;
            add_atom(-1, medium.get_point(node), medium.get_coordination(node));
        }

        elfield.push_back(Vec3(ef[0], ef[1], ef[2]));
        elfield_norm.push_back(ef.norm());
        potential.push_back(pot);
    }
}

const vector<int> SolutionReader::get_node2face_map(Mesh &mesh, int node) {
    expect(node >= 0 && node < mesh.n_nodes_per_face, "Invalid node index!");

    const int n_nodes = mesh.get_n_nodes();
    const int n_faces = mesh.get_n_faces();
    vector<int> map(n_nodes, -1);
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace face = mesh.get_simpleface(i);
        map[face[node]] = i;
    }
    return map;
}

const vector<int> SolutionReader::get_node2elem_map(Mesh &mesh, int node) {
    expect(node >= 0 && node < mesh.n_nodes_per_elem, "Invalid node index!");

    const int n_nodes = mesh.get_n_nodes();
    const int n_elems = mesh.get_n_elems();
    vector<int> map(n_nodes, -1);
    for (int i = 0; i < n_elems; ++i) {
        SimpleElement elem = mesh.get_simpleelem(i);
        map[elem[node]] = i;
    }
    return map;
}

const void SolutionReader::extract_statistics(Mesh &mesh) {
    int i, n_quality, n_nodes;

    vector<int> node2face1 = get_node2face_map(mesh, 0);
    vector<int> node2face2 = get_node2face_map(mesh, 1);
    vector<int> node2face3 = get_node2face_map(mesh, 2);

    vector<int> node2elem1 = get_node2elem_map(mesh, 0);
    vector<int> node2elem2 = get_node2elem_map(mesh, 1);
    vector<int> node2elem3 = get_node2elem_map(mesh, 2);
    vector<int> node2elem4 = get_node2elem_map(mesh, 3);

    n_nodes = mesh.get_n_nodes();
    face_qualities.reserve(n_nodes);
    elem_qualities.reserve(n_nodes);

    mesh.calc_qualities_byface();
    for (i = 0; i < n_nodes; ++i) {
        double q1, q2, q3;
        if ((node2face1[i] >= 0) && (node2face2[i] >= 0) && (node2face3[i] >= 0)) {
            q1 = mesh.get_quality(node2face1[i]);
            q2 = mesh.get_quality(node2face2[i]);
            q3 = mesh.get_quality(node2face3[i]);
        } else {
            q1 = q2 = q3 = 0;
        }
        face_qualities.push_back(max(q1, max(q2, q3)));
    }

    mesh.calc_qualities_byelem();
    for (i = 0; i < n_nodes; ++i) {
        double q1, q2, q3, q4;
        if((node2elem1[i] >= 0) && (node2elem2[i] >= 0) && (node2elem3[i] >= 0) && (node2elem4[i] >= 0)) {
            q1 = mesh.get_quality(node2elem1[i]);
            q2 = mesh.get_quality(node2elem2[i]);
            q3 = mesh.get_quality(node2elem3[i]);
            q4 = mesh.get_quality(node2elem4[i]);
        } else {
            q1 = q2 = q3 = q4 = 0;
        }
        elem_qualities.push_back(max(q1, max(q2, max(q3, q4))));
    }
}

const void SolutionReader::export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez,
        double* Enorm) {
    int i;
    expect(n_atoms >= solution.id.size(), "Solution vector longer than requested!")

    // Initially pass the non-amplified electric field for all the atoms
    for (i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = longrange_efield;
        Enorm[i] = fabs(longrange_efield);
    }

    // Pass the the real electric field for atoms that were in the mesh
    for (i = 0; i < get_n_atoms(); ++i) {
        int identifier = id[i];
        if (identifier < 0 || identifier >= n_atoms) continue;

        Ex[identifier] = elfield[i].x;
        Ey[identifier] = elfield[i].y;
        Ez[identifier] = elfield[i].z;
        Enorm[identifier] = elfield_norm[i];
    }
}
// Reserve memory for solution vectors
const void SolutionReader::reserve(const int n_nodes) {
    Medium::reserve(n_nodes);
    elfield.reserve(n_nodes);
    elfield_norm.reserve(n_nodes);
    potential.reserve(n_nodes);
}

// Map the indices of nodes to the indices of the elements
const vector<int> SolutionReader::get_node2elem_map() {
    vector<int> map(fem->get_n_nodes());
    typename DoFHandler<DIM>::active_cell_iterator cell;

    // Loop through all the elements
    for (cell = fem->dof_handler.begin_active(); cell != fem->dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (unsigned int vertex = 0; vertex < fem->n_verts_per_elem; ++vertex)
            map[cell->vertex_index(vertex)] = cell->active_cell_index();

    return map;
}

// Map the indices of nodes to the indices of the vertex of elements
const vector<int> SolutionReader::get_node2vert_map() {
    vector<int> map(fem->get_n_nodes());
    typename DoFHandler<DIM>::active_cell_iterator cell;

    // Loop through all the elements
    for (cell = fem->dof_handler.begin_active(); cell != fem->dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (unsigned int vertex = 0; vertex < fem->n_verts_per_elem; ++vertex)
            map[cell->vertex_index(vertex)] = vertex;

    return map;
}

// Get mapping between Medium atoms and DealII mesh nodes
const vector<int> SolutionReader::get_medium2node_map(Medium &medium) {
    int i;
    int n_atoms = medium.get_n_atoms();
    vector<int> map(n_atoms, -1);

    typename Triangulation<DIM>::active_vertex_iterator vert;
    vert = fem->triangulation.begin_active_vertex();

    // Loop through mesh vertices that potentially could be on the surface
    for (int vert_cntr = 0; vert_cntr < n_atoms; ++vert_cntr, ++vert)
        // Loop through surface atoms
        for (i = 0; i < n_atoms; ++i) {
            // If map entry already available, take next atom
            if (map[i] >= 0) continue;
            if (medium.get_point(i) == vert->vertex()) map[i] = vert->vertex_index();
        }

    return map;
}

const void SolutionReader::output(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string ftype = get_file_type(file_name);
    expect(ftype == "xyz", "Unsupported file type!");

    int n_atoms = get_n_atoms();

    ofstream out_file(file_name);
    expect(out_file.is_open(), "Can't open a file " + file_name);

    out_file << n_atoms << "\n";
    out_file << get_data_string(-1) << endl;

    for (int i = 0; i < n_atoms; ++i)
        out_file << get_data_string(i) << endl;

    out_file.close();
}

// Compile data string from the data vectors
const string SolutionReader::get_data_string(const int i) {
    if (i < 0)
        return "Solution of DealII operations: id x y z coordination Ex Ey Ez Enorm potential face_quality elem_quality";

    ostringstream strs;
    strs << id[i] << " " << point[i] << " " << coordination[i] << " " << elfield[i] << " "
            << elfield_norm[i] << " " << potential[i] << " " << face_qualities[i] << " "
            << elem_qualities[i];
    return strs.str();
}

} /* namespace femocs */
