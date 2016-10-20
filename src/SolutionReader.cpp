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
SolutionReader::SolutionReader() : longrange_elfield(0), mesh(NULL) {
    reserve(0);
}

SolutionReader::SolutionReader(Mesh* mesh) : longrange_elfield(0), mesh(mesh) {
    reserve(0);
}

const Solution SolutionReader::get_solution(const int i) const {
    expect(i >= 0 && i < get_n_atoms(), "Index out of bounds: " + to_string(i));
    return solution[i];
}

// Return the mapping between atoms & nodes, nodes & elements and nodes & vertices
const void SolutionReader::get_maps(DealII& fem, vector<int>& medium2node, vector<int>& node2elem, vector<int>& node2vert) {
    const int n_atoms = mesh->get_n_nodes();
    const int n_nodes = fem.get_n_nodes();

    medium2node.resize(n_atoms, -1);
    node2elem.resize(n_nodes);
    node2vert.resize(n_nodes);

    typename Triangulation<DIM>::active_vertex_iterator vertex;
    typename DoFHandler<DIM>::active_cell_iterator cell;
    int j;

    // Loop through mesh vertices that potentially could be on the Medium
    for (j = 0, vertex = fem.triangulation.begin_active_vertex(); j < n_atoms; ++j, ++vertex)
        // Loop through Medium atoms
        for (int i = 0; i < n_atoms; ++i)
            if ( (medium2node[i] < 0) && (mesh->get_node(i) == vertex->vertex()) ) {
                medium2node[i] = vertex->vertex_index();
                break;
            }

    // Loop through all the mesh elements
    for (cell = fem.dof_handler.begin_active(); cell != fem.dof_handler.end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < fem.n_verts_per_elem; ++i) {
            node2elem[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }
}

// Extract the electric potential and electric field values on Medium atoms from FEM solution
const void SolutionReader::extract_solution(DealII &fem) {
    longrange_elfield = fem.get_elfield();
    const int n_nodes = mesh->get_n_nodes();
    reserve(n_nodes);

    vector<int> medium2node, node2elem, node2vert;
    get_maps(fem, medium2node, node2elem, node2vert);

    vector<int> cell_indxs; cell_indxs.reserve(n_nodes);
    vector<int> vert_indxs; vert_indxs.reserve(n_nodes);

    for (int n : medium2node)
        if (n >= 0) {
            cell_indxs.push_back(node2elem[n]);
            vert_indxs.push_back(node2vert[n]);
        }

    vector<Vec3> ef = fem.get_elfield(cell_indxs, vert_indxs);
    vector<double> pot = fem.get_potential(cell_indxs, vert_indxs);

    int i = 0;
    for (int node = 0; node < n_nodes; ++node) {
        add_atom( mesh->get_node(node) );

        if (medium2node[node] >= 0) {
            solution.push_back( Solution(ef[i], ef[i].length(), pot[i]) );
            i++;
        }
        else
            solution.push_back( Solution(error_field) );
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
        skip_atom[i] = solution[i].el_norm <= fabs(longrange_elfield);

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
        skip_atom[i] = (solution[i].el_norm <= fabs(longrange_elfield)) || (solution[i].el_norm == error_field);

    // Average by moving up
    for (int i = 0; i < n_atoms-1; ++i) {
        int i_active = i;
        int i_neighb = i+1;

        if (skip_atom[i_active] || skip_atom[i_neighb]) continue;

        solution[i_active].elfield = get_ema(i_active, i_neighb, smooth_width);
        solution[i_active].el_norm = solution[i_active].elfield.length();
    }
}

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

    // Sort atoms back to their initial order
    sort_atom_id();
}

const void SolutionReader::get_histogram(vector<int> &bins, vector<double> &bounds, const int coord) {
//    require(coord >= 0 && coord <= 3, "Invalid component: " + to_string(coord));

    const int n_atoms = get_n_atoms();
    const int n_bins = bins.size();
    const int n_bounds = bounds.size();
    const double eps = 1e-5;

    // Find minimum and maximum values from all non-error values
    double value_min = 1e100;
    double value_max =-1e100;
    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coord == 3) value = abs(solution[i].el_norm);
        else            value = abs(solution[i].elfield[coord]);

        if (value < error_field) {
            value_min = min(value_min, value);
            value_max = max(value_max, value);
        }
    }

    // Fill the bounds with values value_min:value_step:(value_max + epsilon)
    // Epsilon is added to value_max to include the maximum value in the up-most bin
    double value_step = (value_max - value_min) / n_bins;
    for (int i = 0; i < n_bounds; ++i)
        bounds[i] = value_min + value_step * i;
    bounds[n_bounds-1] += eps;

    for (int i = 0; i < n_atoms; ++i)
        for (int j = 0; j < n_bins; ++j) {
            if (coord == 3) value = abs(solution[i].el_norm);
            else            value = abs(solution[i].elfield[coord]);

            if (value >= bounds[j] && value < bounds[j+1]) {
                bins[j]++;
                continue;
            }
        }
}

// Return the element index that contains n-th node
const int SolutionReader::get_elem(const int n) {
    for (int i = 0; i < mesh->get_n_elems(); ++i)
        if (mesh->get_simpleelem(i) == n)
            return i;
    return -1;
}

const Vec3 SolutionReader::get_average_solution(const int I) {
    int elem = get_elem(I);                  // get some element that contains I-th node
    if (elem < 0) return Vec3(error_field);  // if no such element exists, return error field

    vector<SimpleElement> selems;
    selems.push_back(mesh->get_simpleelem(elem));
    for (int nbr : mesh->get_elem_neighbours(elem))
        if (nbr >= 0 && mesh->get_simpleelem(nbr) == I)
            selems.push_back(mesh->get_simpleelem(nbr));

    Vec3 average(0.0);
    int n_averages = 0;

    // Average the node value over other nodes of the element
    for (SimpleElement selem : selems)
        for (int i = 0; i < mesh->n_nodes_per_elem; ++i)
            if (selem[i] != I) {
                average += solution[selem[i]].elfield;
                n_averages++;
            }

    average *= (1.0 / n_averages);
    return average;
}

// Function to clean the result from peaks
const void SolutionReader::clean(const int coordinate, const int n_bins) {
//    require(coordinate >= 0 && coordinate <= 3, "Invalid coordinate: " + to_string(coordinate));
    if (n_bins <= 1) return;

    const int n_atoms = get_n_atoms();

    vector<int> bins(n_bins, 0);
    vector<double> bounds(n_bins+1);
    get_histogram(bins, bounds, coordinate);

    // Find the first bin with zero entries; this will determine the maximum allowed elfield norm
    double value_max = bounds[n_bins];
    for (int i = n_bins-1; i >= 0; --i)
        if (bins[i] == 0) value_max = bounds[i];

    cout.precision(3);
    cout << endl << coordinate << " " << value_max << endl;
    for (int i = 0; i < bins.size(); ++i)
        cout << bins[i] << " ";
    cout << endl;
    for (int i = 0; i < bounds.size(); ++i)
        cout << bounds[i] << " ";
    cout << endl;

    // If all the bins are filled, no blocking will be applied
    if (value_max == bounds[n_bins])
        return;

    double value;
    for (int i = 0; i < n_atoms; ++i) {
        if (coordinate == 3) value = abs(solution[i].el_norm);
        else                 value = abs(solution[i].elfield[coordinate]);

        if (value > value_max) {
            solution[i].elfield = get_average_solution(i);
            solution[i].el_norm = solution[i].elfield.length();
        }
    }
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
