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
SolutionReader::SolutionReader() : mesh(NULL) {
    reserve(0);
}

SolutionReader::SolutionReader(TetgenMesh* mesh) : mesh(mesh) {
    reserve(0);
}

// Return pointer to the tetrahedral mesh SolutionReader is using
TetgenMesh* SolutionReader::get_mesh() {
    return mesh;
}

const Solution SolutionReader::get_solution(const int i) const {
    expect(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return solution[i];
}

/* Return the mapping between tetrahedral & hexahedral mesh nodes, 
   nodes & hexahedral elements and nodes & element's vertices  */
const void SolutionReader::get_maps(DealII& fem, vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert) {
    const int n_tet_nodes = mesh->nodes.size();
    const int n_hex_nodes = fem.triangulation.n_used_vertices();

    tet2hex.resize(n_tet_nodes, -1);
    node2hex.resize(n_hex_nodes);
    node2vert.resize(n_hex_nodes);

    // Loop through hexahedral mesh vertices
    typename Triangulation<DIM>::active_vertex_iterator vertex = fem.triangulation.begin_active_vertex();
    for (int j = 0; j < n_tet_nodes; ++j, ++vertex)
        // Loop through tetrahedral mesh vertices
        for (int i = 0; i < n_tet_nodes; ++i)
            if ( (tet2hex[i] < 0) && (mesh->nodes[i] == vertex->vertex()) ) {
                tet2hex[i] = vertex->vertex_index();
                break;
            }

    // Loop through all the hexahedral mesh elements
    typename DoFHandler<DIM>::active_cell_iterator cell;
    for (cell = fem.dof_handler.begin_active(); cell != fem.dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < fem.n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
const void SolutionReader::extract_solution(DealII &fem) {
    const int n_nodes = mesh->nodes.size();
    reserve(n_nodes);

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tet2hex, node2hex, node2vert;
    get_maps(fem, tet2hex, node2hex, node2vert);

    vector<int> cell_indxs; cell_indxs.reserve(n_nodes);
    vector<int> vert_indxs; vert_indxs.reserve(n_nodes);

    for (int n : tet2hex)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }

    vector<Vec3> ef = fem.get_elfield(cell_indxs, vert_indxs);      // electric field
    vector<double> pot = fem.get_potential(cell_indxs, vert_indxs); // electric potential

    int i = 0;
    for (int node = 0; node < n_nodes; ++node) {
        add_atom( mesh->nodes[node] );

        // If there is a common node between tet and hex meshes, store actual solution
        if (tet2hex[node] >= 0) {
            solution.push_back( Solution(ef[i], ef[i].norm(), pot[i]) );
            i++;
        }
        
        // In case of non-common node, store solution with error value
        else
            solution.push_back( Solution(error_field) );
    }
}

// Sort solution vectors by atom id
const void SolutionReader::sort_atom_id() {
    const int n_atoms = get_n_atoms();

    // Sort atoms by their ID
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].sort_indx = i;
    sort(atoms.begin(), atoms.end(), Atom::sort_id());

    // Sort solution vectors into same order as atoms
    for (int i = 0; i < n_atoms; ++i)
        solution[atoms[i].sort_indx].sort_indx = i;
    sort( solution.begin(), solution.end(), Solution::sort_up() );
}

// Sort the atoms and results
const void SolutionReader::sort_atoms(const int x1, const int x2, const string& direction) {
    const int n_atoms = get_n_atoms();

    // Sort atoms into desired order
    for (int i = 0; i < n_atoms; ++i)
        atoms[i].sort_indx = i;
    Medium::sort_atoms(x1, x2, direction);

    // Sort solution vectors into same order as atoms
    for (int i = 0; i < n_atoms; ++i)
        solution[atoms[i].sort_indx].sort_indx = i;

    if (direction == "up" || direction == "asc")
        sort( solution.begin(), solution.end(), Solution::sort_up() );
    else if (direction == "down" || direction == "desc")
        sort( solution.begin(), solution.end(), Solution::sort_down() );
}

// Print some statistics about the top region of tip
const void SolutionReader::print_statistics() {
#if not VERBOSE
    return;
#endif

    calc_statistics();

    double zmax = sizes.zmax - 3.0;
    double mn = 1e20; double mx = -1e20; double avg = 0; int cntr = 0;

    for (int i = 0; i < get_n_atoms(); ++i)
        if (get_point(i).z > zmax && solution[i].el_norm < error_field) {
            mn = min(mn, solution[i].el_norm);
            mx = max(mx, solution[i].el_norm);
            avg += solution[i].el_norm;
            cntr++;
        }

    cout << "zmax, n_atoms:\t" << zmax << ", " << cntr
            << "\nmin, max, span:\t" << mn << ", " << mx << ", " << mx-mn <<
            "\naverage:\t" << avg / cntr << endl;
}

// Reserve memory for solution vectors
const void SolutionReader::reserve(const int n_nodes) {
    atoms.clear();
    solution.clear();

    Medium::reserve(n_nodes);
    solution.reserve(n_nodes);
}

// Compile data string from the data vectors for file output
const string SolutionReader::get_data_string(const int i) {
    if (i < 0) return "SolutionReader data: id x y z coordination Ex Ey Ez Enorm potential";

    ostringstream strs;
//    strs << i << " " << atoms[i].point << " 0" << " " << solution[i];
    strs << atoms[i] << " " << solution[i];
    return strs.str();
}

} /* namespace femocs */
