/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include "Macros.h"

#include <fstream>
#include <float.h>

using namespace std;
namespace femocs {

// SolutionReader constructor
SolutionReader::SolutionReader() {
    reserve(0);
}

const Solution SolutionReader::get_solution(const int i) const {
    expect(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return solution[i];
}

/* Return the mapping between tetrahedral & hexahedral mesh nodes, 
   nodes & hexahedral elements and nodes & element's vertices  */
const void SolutionReader::get_maps(DealII& fem, vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert) {
    const int n_tet_nodes = get_n_atoms();
    const int n_hex_nodes = fem.triangulation.n_used_vertices();

    tet2hex.resize(n_tet_nodes, -1);
    node2hex.resize(n_hex_nodes);
    node2vert.resize(n_hex_nodes);

    // Loop through the hexahedral mesh vertices
    typename Triangulation<DIM>::active_vertex_iterator vertex = fem.triangulation.begin_active_vertex();
    for (int j = 0; j < n_tet_nodes; ++j, ++vertex)
        // Loop through tetrahedral mesh vertices
        for (int i = 0; i < n_tet_nodes; ++i)
            if ( (tet2hex[i] < 0) && (get_point(i) == vertex->vertex()) ) {
                tet2hex[i] = vertex->vertex_index();
                break;
            }

    // Loop through the hexahedral mesh elements
    typename DoFHandler<DIM>::active_cell_iterator cell;
    for (cell = fem.dof_handler.begin_active(); cell != fem.dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < fem.n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
const void SolutionReader::extract_solution(DealII &fem, const TetgenNodes &nodes) {
    const int n_nodes = nodes.size();

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int node = 0; node < n_nodes; ++node)
        add_atom( Atom(node, nodes[node], -1) );

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tet2hex, node2hex, node2vert;
    get_maps(fem, tet2hex, node2hex, node2vert);

    // Generate lists of hexahedra and hexahedra nodes where the tetrahedra nodes are located
    vector<int> cell_indxs; cell_indxs.reserve(n_nodes);
    vector<int> vert_indxs; vert_indxs.reserve(n_nodes);
    for (int n : tet2hex)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }

    vector<Vec3> ef = fem.get_elfield(cell_indxs, vert_indxs);      // get list of electric fields
    vector<double> pot = fem.get_potential(cell_indxs, vert_indxs); // get list of electric potentials

    require( ef.size() == pot.size(), "Mismatch of vector sizes: "
            + to_string(ef.size())  + ", " + to_string(pot.size()) );

    int i = 0;
    for (int node = 0; node < n_nodes; ++node) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (tet2hex[node] >= 0) {
            require(i < ef.size(), "Invalid index: " + to_string(i));
            solution.push_back( Solution(ef[i], ef[i].norm(), pot[i]) );
            i++;
        }
        
        // In case of non-common node, store solution with error value
        else
            solution.push_back( Solution(error_field) );
    }
}

// Reserve memory for solution vectors
const void SolutionReader::reserve(const int n_nodes) {
    require(n_nodes >= 0, "Invalid number of nodes: " + to_string(n_nodes));

    atoms.clear();
    solution.clear();

    atoms.reserve(n_nodes);
    solution.reserve(n_nodes);
}

// Compile data string from the data vectors for file output
const string SolutionReader::get_data_string(const int i) {
    if (i < 0) return "SolutionReader data: id x y z dummy Ex Ey Ez Enorm potential";

    ostringstream strs;
    strs << atoms[i] << " " << solution[i];
    return strs.str();
}

} // namespace femocs
