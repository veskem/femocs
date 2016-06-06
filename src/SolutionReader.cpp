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
SolutionReader::SolutionReader() : longrange_efield(0) {
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

        elfield.push_back( Vec3(ef[0], ef[1], ef[2]) );
        elfield_norm.push_back(ef.norm());
        potential.push_back(pot);
    }
}

const void SolutionReader::export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
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
        int ident = id[i];
        if(ident < 0 || ident >= n_atoms) continue;

        Ex[ident] = elfield[i].x;
        Ey[ident] = elfield[i].y;
        Ez[ident] = elfield[i].z;
        Enorm[ident] = elfield_norm[i];
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

// Compile data string from the data vectors
const string SolutionReader::get_data_string(const int i) {
    if(i < 0)
        return "Solution of DealII operations: id x y z coordination Ex Ey Ez Enorm potential";

    ostringstream strs;
    strs << id[i] << " " << point[i] << " " << coordination[i]
            << " " << elfield[i] << " " << elfield_norm[i] << potential[i];
    return strs.str();
}

} /* namespace femocs */
