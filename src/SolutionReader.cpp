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
SolutionReader::SolutionReader() :
        longrange_efield(0) {
    reserve(0);
}

const void SolutionReader::extract_solution(DealII* fem, Medium &medium) {
    dealii::Tensor<1, DIM> ef;
    double pot;

    double t0, t1, tstart;

    tstart = omp_get_wtime();

    this->fem = fem;
    longrange_efield = fem->get_elfield();
    int n_nodes = medium.get_n_atoms();
    reserve(n_nodes);

    vector<int> medium2node_map = get_medium2node_map(medium);
    vector<int> node2elem_map = get_node2elem_map();
    vector<int> node2vert_map = get_node2vert_map();

    t0 = omp_get_wtime() - tstart;

    for (int node = 0; node < n_nodes; ++node) {
        int n = medium2node_map[node];
        if (n >= 0) {
            ef = fem->get_elfield_at_node(node2elem_map[n], node2vert_map[n]);
            pot = fem->get_potential_at_node(node2elem_map[n], node2vert_map[n]);
            add_atom(medium.get_id(node), medium.get_point(node), medium.get_coordination(node));
        } else {
            ef[0] = ef[1] = ef[2] = pot = error_field;
            add_atom(-1, medium.get_point(node), medium.get_coordination(node));
        }

        elfield.push_back(Vec3(ef[0], ef[1], ef[2]));
        elfield_norm.push_back(ef.norm());
        potential.push_back(pot);
    }

    t1 = omp_get_wtime() - t0 - tstart;
    tstart = omp_get_wtime() - tstart;

    cout << "\npercentage 0: " << t0/tstart << "\npercentage 1: " << t1/tstart << endl;
}

const void SolutionReader::extract_solution_vol2(DealII* fem, Medium &medium) {
    double t0, t1, tstart;
    tstart = omp_get_wtime();

    this->fem = fem;
    longrange_efield = fem->get_elfield();
    int n_nodes = medium.get_n_atoms();
    reserve(n_nodes);

    vector<int> medium2node_map = get_medium2node_map(medium);
    vector<int> node2elem_map = get_node2elem_map();
    vector<int> node2vert_map = get_node2vert_map();

    t0 = omp_get_wtime() - tstart;


    elfield.resize(n_nodes, Vec3(error_field));
    potential.resize(n_nodes, error_field);
    elfield_norm.resize(n_nodes, error_field);

    vector<int> cell_indxs;
    vector<int> vert_indxs;
    cell_indxs.reserve(n_nodes);
    vert_indxs.reserve(n_nodes);

    for (int n : medium2node_map)
        if (n >= 0) {
            cell_indxs.push_back(node2elem_map[n]);
            vert_indxs.push_back(node2vert_map[n]);
        }

    vector<Vec3> ef = fem->get_elfield_at_node(cell_indxs, vert_indxs);
//    vector<double> pot = fem->get_potential_at_node(cell_indxs, vert_indxs);

    int indx = 0;

    for (int node = 0; node < n_nodes; ++node)
        if (medium2node_map[node] >= 0) {
//            potential[node] = pot[indx];
            elfield[node] = ef[indx++];
            elfield_norm[node] = elfield[node].length();
            add_atom(medium.get_id(node), medium.get_point(node), medium.get_coordination(node));
        } else {
            add_atom(-1, medium.get_point(node), medium.get_coordination(node));
        }

    t1 = omp_get_wtime() - t0 - tstart;
    tstart = omp_get_wtime() - tstart;

    cout << "\npercentage 0: " << t0/tstart << "\npercentage 1: " << t1/tstart << endl;
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

inline Vec3 SolutionReader::get_sma_up(const int i, const int n_samples) {
    Vec3 movavg(0.0);
    for(int j = i; j < i + n_samples ; ++j)
        movavg += elfield[j];

    return movavg / n_samples;
}

inline Vec3 SolutionReader::get_sma_down(const int i, const int n_samples) {
    Vec3 movavg(0.0);
    for(int j = i; j > i - n_samples ; --j)
        movavg += elfield[j];

    return movavg / n_samples;
}

inline Vec3 SolutionReader::get_sma_up(const int i, const int n_samples, const vector<int>& map) {
    Vec3 movavg(0.0);
    for(int j = i; j < i + n_samples ; ++j)
        movavg += elfield[map[j]];

    return movavg / n_samples;
}

inline Vec3 SolutionReader::get_sma_down(const int i, const int n_samples, const vector<int>& map) {
    Vec3 movavg(0.0);
    for(int j = i; j > i - n_samples ; --j)
        movavg += elfield[map[j]];

    return movavg / n_samples;
}

inline Vec3 SolutionReader::get_ema(const int i_active, const int i_neighb, const int n_average) {
    const double k = 2.0 / (n_average + 1);
    return (elfield[i_active] * k) + (elfield[i_neighb] * (1 - k));
}

// Function to get sorting indices for atoms that are sorted by their radial coordinates
const vector<int> SolutionReader::get_radial_map() {
    const int n_atoms = get_n_atoms();

    calc_statistics();

    // Calculate the vector with atom distances from the centre axis of simubox
    vector<double> radiuses(n_atoms);
    Point2 origin((sizes.xmax + sizes.xmin) / 2, (sizes.ymax + sizes.ymin) / 2);
    for(int i = 0; i < n_atoms; ++i)
        radiuses[i] = origin.distance2(get_point2(i));

    // Generate vector with indices [0, n_atoms-1]
    vector<int> indxs(n_atoms);
    size_t n(0);
    generate(indxs.begin(), indxs.end(), [&]{ return n++; });

    // Sort indexes by the radiuses
    auto comparator = [&radiuses](int a, int b){ return radiuses[a] < radiuses[b]; };
    sort(indxs.begin(), indxs.end(), comparator);

    return indxs;
}

// Function to get sorting indices for atoms that are sorted by their x- and y-coordinates
const vector<int> SolutionReader::get_cartesian_map(const int x1, const int x2) {
    expect((x1 == 0 || x1 == 1) && (x2 == -1 || x2 == 0 || x2 == 1), "Invalid arguments!");

    // Generate vector with indices [0, n_atoms-1]
    vector<int> indxs(get_n_atoms());
    size_t n(0);
    generate(indxs.begin(), indxs.end(), [&]{ return n++; });

    // Generate duplicate of points
    vector<Point3> points = point;

    // Sort indices by x-coordinates
    auto comparator_x = [&points](int a, int b) { return points[a].x < points[b].x; };
    // Sort indices by y-coordinates
    auto comparator_y = [&points](int a, int b) { return points[a].y < points[b].y; };

    // Sort indices first by x and then by y coordinates
    auto comparator_xy = [&points](int a, int b) {
        return (points[a].x < points[b].x) || ((points[a].x == points[b].x) && (points[a].y < points[b].y)); };
    // Sort indices first by y and then by x coordinates
    auto comparator_yx = [&points](int a, int b) {
        return (points[a].y < points[b].y) || ((points[a].y == points[b].y) && (points[a].x < points[b].x)); };

    // Sort indices according to function arguments
    if (x1 == 0 && x2 == 1)
        sort(indxs.begin(), indxs.end(), comparator_xy);
    else if (x1 == 1 && x2 == 0)
        sort(indxs.begin(), indxs.end(), comparator_yx);
    else if (x1 == 0 && (x2 == 0 || x2 == -1))
        sort(indxs.begin(), indxs.end(), comparator_x);
    else if (x1 == 1 && (x2 == 1 || x2 == -1))
        sort(indxs.begin(), indxs.end(), comparator_y);

    return indxs;
}

const void SolutionReader::smoothen_result(const int n_average, const int repetitions) {
    if (n_average <= 0 || repetitions <= 0) return;

    // Repeat the averaging cycles for desired times to get properly centred results
    for (int r = 0; r < repetitions; ++r) {
        smoothen_result_ema(n_average, get_cartesian_map(1, 0));
        smoothen_result_ema(n_average, get_cartesian_map(0, 1));
    }

#if VERBOSE
    // Print some statistics about the top region of tip
    calc_statistics();

    double zmax = sizes.zmax - 3.0;// 33.0;
    double mn = 1e20; double mx = -1e20; double avg = 0; int cntr = 0;

    for (int i = 0; i < get_n_atoms(); ++i)
        if (point[i].z > zmax) {
            mn = min(mn, elfield_norm[i]);
            mx = max(mx, elfield_norm[i]);
            avg += elfield_norm[i];
            cntr++;
        }

    cout << "\n\nzmax, n_atoms, n_average: " << zmax << ", " << cntr << ", " << n_average
            << "\nmin, max: " << mn << ", " << mx
            << "\nspan: " << mx-mn << "\naverage: " << avg / cntr << "\n\n";
#endif
}

const void SolutionReader::smoothen_result_sma(const int n_average) {
    // Atoms from non-enhanced electric field regions will be skipped from averaging
    const vector<bool> skip_atom = vector_less_equal(&elfield_norm, fabs(longrange_efield));
    const int n_atoms = get_n_atoms();

    // Average by moving up
    // moving up is more sensitive, thus the division by two of n_average
    for (int i = 0; i < (n_atoms - n_average); ++i) {
        if (skip_atom[i]) continue;

        elfield[i] = get_sma_up(i, n_average);
        elfield_norm[i] = elfield[i].length();
    }
    // Average by moving down
    for (int i = (n_atoms - 1); i >= (n_average - 1); --i) {
        if (skip_atom[i]) continue;

        elfield[i] = get_sma_down(i, n_average);
        elfield_norm[i] = elfield[i].length();
    }
}

const void SolutionReader::smoothen_result_sma(const int n_average, const vector<int>& sort_map) {
    // Atoms from non-enhanced electric field regions will be skipped from averaging
    const vector<bool> skip_atom = vector_less_equal(&elfield_norm, fabs(longrange_efield));
    const int n_atoms = get_n_atoms();

    // Average by moving up
    // moving up is more sensitive, thus the division by two of n_average
    for (int i = 0; i < (n_atoms - n_average); ++i) {
        int I = sort_map[i];
        if (skip_atom[I]) continue;
        elfield[I] = get_sma_up(i, n_average, sort_map);
        elfield_norm[I] = elfield[I].length();
    }
//    // Average by moving down
//    for (int i = (n_atoms - 1); i >= (n_average - 1); --i) {
//        int I = sort_map[i];
//        if (skip_atom[I]) continue;
//        elfield[I] = get_sma_down(i, n_average, sort_map);
//        elfield_norm[I] = elfield[I].length();
//    }
}

const void SolutionReader::smoothen_result_ema(const int n_average) {
    // Atoms from non-enhanced electric field regions will be skipped from averaging
    const vector<bool> skip_atom = vector_less_equal(&elfield_norm, fabs(longrange_efield));
    const int n_atoms = get_n_atoms();

    // Average by moving up
    for (int i = 1; i < n_atoms; ++i) {
        int i_active = i;
        int i_neighb = i-1;

        if (skip_atom[i_active]) continue;

        elfield[i_active] = get_ema(i_active, i_neighb, n_average);
        elfield_norm[i_active] = elfield[i_active].length();
    }
    // Average by moving down
//    for (int i = (n_atoms - 2); i >= 0; --i) {
//        int i_active = i;
//        int i_neighb = i+1;
//
//        if (skip_atom[i_active]) continue;
//
//        elfield[i_active] = get_ema(i_active, i_neighb, n_average);
//        elfield_norm[i_active] = elfield[i_active].length();
//    }
}

const void SolutionReader::smoothen_result_ema(const int n_average, const vector<int>& sort_map) {
    // Atoms from non-enhanced electric field regions will be skipped from averaging
    const vector<bool> skip_atom = vector_less_equal(&elfield_norm, fabs(longrange_efield));
    const int n_atoms = get_n_atoms();

    // Average by moving up
    for (int i = 1; i < n_atoms; ++i) {
        int i_active = sort_map[i];
        int i_neighb = sort_map[i-1];

        if (skip_atom[i_active]) continue;

        elfield[i_active] = get_ema(i_active, i_neighb, n_average);
        elfield_norm[i_active] = elfield[i_active].length();
    }
    // Average by moving down
//    for (int i = (n_atoms - 2); i >= 0; --i) {
//        int i_active = sort_map[i];
//        int i_neighb = sort_map[i+1];
//
//        if (skip_atom[i_active]) continue;
//
//        elfield[i_active] = get_ema(i_active, i_neighb, n_average);
//        elfield_norm[i_active] = elfield[i_active].length();
//    }
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

// Compile data string from the data vectors
const string SolutionReader::get_data_string(const int i) {
    if (i < 0) return "SolutionReader data: id x y z coordination Ex Ey Ez Enorm potential";// face_quality elem_quality";

    ostringstream strs;
//    strs << id[i] << " " << point[i] << " " << coordination[i] << " " << elfield[i] << " "
    strs << i << " " << point[i] << " " << id[i] << " " << elfield[i] << " "
            << elfield_norm[i] << " " << potential[i];
            // << " " << face_qualities[i] << " " << elem_qualities[i];
    return strs.str();
}

} /* namespace femocs */
