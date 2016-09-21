/*
 * SolutionReader.cpp
 *
 *  Created on: 6.6.2016
 *      Author: veske
 */

#include "SolutionReader.h"
#include <fstream>

#include <omp.h>

using namespace std;
namespace femocs {

// SolutionReader constructor
SolutionReader::SolutionReader() : longrange_efield(0) {
    reserve(0);
}

// Return the mapping between atoms & nodes, nodes & elements and nodes & vertices
const void SolutionReader::get_maps(Medium& medium, vector<int>& medium2node, vector<int>& node2elem, vector<int>& node2vert) {
    const int n_atoms = medium.get_n_atoms();
    const int n_nodes = fem->get_n_nodes();

    medium2node.resize(n_atoms, -1);
    node2elem.resize(n_nodes);
    node2vert.resize(n_nodes);

    typename Triangulation<DIM>::active_vertex_iterator vertex;
    typename DoFHandler<DIM>::active_cell_iterator cell;

    int i, j;

    // Loop through mesh vertices that potentially could be on the Medium
    for (j = 0, vertex = fem->triangulation.begin_active_vertex(); j < n_atoms; ++j, ++vertex)
        // Loop through Medium atoms
        for (i = 0; i < n_atoms; ++i)
            if ( (medium2node[i] < 0) && (medium.get_point(i) == vertex->vertex()) ) {
                medium2node[i] = vertex->vertex_index();
                break;
            }

    // Loop through all the mesh elements
    for (cell = fem->dof_handler.begin_active(); cell != fem->dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (i = 0; i < fem->n_verts_per_elem; ++i) {
            node2elem[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
}

// Extract the electric potential and electric field values on Medium atoms from FEM solution
const void SolutionReader::extract_solution(DealII* fem, Medium &medium) {
    this->fem = fem;
    longrange_efield = fem->get_elfield();
    const int n_nodes = medium.get_n_atoms();
    reserve(n_nodes);

    vector<int> medium2node, node2elem, node2vert;
    get_maps(medium, medium2node, node2elem, node2vert);

    vector<int> cell_indxs; cell_indxs.reserve(n_nodes);
    vector<int> vert_indxs; vert_indxs.reserve(n_nodes);

    for (int n : medium2node)
        if (n >= 0) {
            cell_indxs.push_back(node2elem[n]);
            vert_indxs.push_back(node2vert[n]);
        }

    vector<Vec3> ef = fem->get_elfield_at_node(cell_indxs, vert_indxs);
    vector<double> pot = fem->get_potential_at_node(cell_indxs, vert_indxs);

    int i, node;
    for (i = 0, node = 0; node < n_nodes; ++node)
        if (medium2node[node] >= 0) {
            solution.push_back( Solution(ef[i], ef[i].length(), pot[i]) );
            add_atom(medium.get_atom(node));
            i++;
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

// Get simple moving average value in up direction
inline Vec3 SolutionReader::get_sma_up(const int i, const int n_samples) {
    Vec3 movavg(0.0);
    for(int j = i; j < i + n_samples ; ++j)
        movavg += solution[j].elfield;

    return movavg / n_samples;
}

// Get simple moving average value in down direction
inline Vec3 SolutionReader::get_sma_down(const int i, const int n_samples) {
    Vec3 movavg(0.0);
    for(int j = i; j > i - n_samples ; --j)
        movavg += solution[j].elfield;

    return movavg / n_samples;
}

// Get exponential moving average value
inline Vec3 SolutionReader::get_ema(const int i_active, const int i_neighb, const double n_average) {
    const double k = 2.0 / (n_average + 1);
    return (solution[i_active].elfield * k) + (solution[i_neighb].elfield * (1 - k));
}

// Smoothen the electric field with simple moving average technique
const void SolutionReader::smoothen_result_sma(const int n_average) {
    const int n_atoms = get_n_atoms();

    // Atoms from non-enhanced electric field regions will be skipped from averaging
    vector<bool> skip_atom(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        skip_atom[i] = solution[i].el_norm <= fabs(longrange_efield);

    // Average by moving up
    // moving up is more sensitive, thus the division by two of n_average
    for (int i = 0; i < (n_atoms - n_average); ++i) {
        if (skip_atom[i]) continue;

        solution[i].elfield = get_sma_up(i, n_average);
        solution[i].el_norm = solution[i].elfield.length();
    }
    // Average by moving down
    for (int i = (n_atoms - 1); i >= (n_average - 1); --i) {
        if (skip_atom[i]) continue;

        solution[i].elfield = get_sma_down(i, n_average);
        solution[i].el_norm = solution[i].elfield.length();
    }
}

// Smoothen the electric field with exponential moving average technique
const void SolutionReader::smoothen_result_ema(const double smooth_width) {
    const int n_atoms = get_n_atoms();

    // Atoms from non-enhanced electric field regions will be skipped from averaging
    vector<bool> skip_atom(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        skip_atom[i] = (solution[i].el_norm <= fabs(longrange_efield)) || (solution[i].el_norm == error_field);

    // Average by moving up
    for (int i = 0; i < n_atoms-1; ++i) {
        int i_active = i;
        int i_neighb = i+1;

        if (skip_atom[i_active] || skip_atom[i_neighb]) continue;

        solution[i_active].elfield = get_ema(i_active, i_neighb, smooth_width);
        solution[i_active].el_norm = solution[i_active].elfield.length();
    }
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

// Function to perform smoothing of electric field
const void SolutionReader::smoothen_result(const double smooth_width) {
    if (smooth_width <= 0) return;

    // Sort atoms first by x- and then by y-coordinate
    sort_atoms(0, 1, "up");
    // Sweep the results in sorted order
    smoothen_result_ema(smooth_width);

    // Sort atoms first by y- and then by x-coordinate
    sort_atoms(1, 0, "up");
    // Sweep the results in sorted order
    smoothen_result_ema(smooth_width);
}

// Print some statistics about the top region of tip
const void SolutionReader::print_statistics() {
#if not VERBOSE
    return;
#endif

    calc_statistics();

    double zmax = sizes.zmax - 3.0;// 33.0;
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

const void SolutionReader::export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez,
        double* Enorm) {
    int i;
    expect(n_atoms >= get_n_atoms(), "Solution vector longer than requested!")

    // Initially pass the long range electric field for all the atoms
    // TODO: THAT FIELD IS ALSO PASSED FOR THE BULK ATOMS (IT SHOULDN'T). FIX THIS!
    for (i = 0; i < n_atoms; ++i) {
        Ex[i] = 0;
        Ey[i] = 0;
        Ez[i] = longrange_efield;
        Enorm[i] = fabs(longrange_efield);
    }

    // Pass the the real electric field for surface atoms that were in the mesh
    for (i = 0; i < get_n_atoms(); ++i) {
        int identifier = get_id(i);
        if (identifier < 0 || identifier >= n_atoms) continue;

        Ex[identifier] = solution[i].elfield.x;
        Ey[identifier] = solution[i].elfield.y;
        Ez[identifier] = solution[i].elfield.z;
        Enorm[identifier] = solution[i].el_norm;
    }
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
    if (i < 0) return "SolutionReader data: id x y z coordination Ex Ey Ez Enorm potential";// face_quality elem_quality";

    ostringstream strs;
//    strs << i << " " << atoms[i].point << " 0" << " " << solution[i];
    strs << atoms[i] << " " << solution[i];
            // << " " << face_qualities[i] << " " << elem_qualities[i];
    return strs.str();
}

} /* namespace femocs */
