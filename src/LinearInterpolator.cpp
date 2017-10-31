/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "Macros.h"
#include "LinearInterpolator.h"

#include "float.h"

using namespace std;
namespace femocs {

/* ==================================================================
 *  ===================== LinearInterpolator =======================
 * ================================================================== */

// Calculate distance-dependent weights for a point with respect to the cell
template<int dim>
array<double,dim> LinearInterpolator<dim>::get_weights(const Point3 &point, const SimpleCell<dim>& scell) const {
    array<double,dim> weights;
    double w_sum = 0;

    for (int i = 0; i < dim; ++i) {
        double w = exp(decay_factor * point.distance(vertices[scell[i]]));
        weights[i] = w;
        w_sum += w;
    }
    // check for zero and nan
    require(w_sum > 0 && w_sum == w_sum, "Invalid interpolation weight: " + to_string(w_sum));
    w_sum = 1.0 / w_sum;

    for (int i = 0; i < dim; ++i)
        weights[i] *= w_sum;

    return weights;
}

// Interpolate both scalar and vector data inside or near the cell
template<int dim>
Solution LinearInterpolator<dim>::interp_solution(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

//    SimpleCell<dim> scell = cells[cell];
    array<unsigned, dim> scell = solution_indices[cell];

    // calculate weights or barycentric coordinates
    array<double,dim> bcc;
    if (c >= 0) bcc = get_bcc(Vec3(point), cell);
    else bcc = get_weights(point, cell);

    // Interpolate electric field
    Vec3 vector_i(0.0);
    for (int i = 0; i < dim; ++i)
        vector_i += solutions[scell[i]].vector * bcc[i];

    // Interpolate potential
    double scalar_i(0.0);
    for (int i = 0; i < dim; ++i)
        scalar_i += solutions[scell[i]].scalar * bcc[i];

    return Solution(vector_i, scalar_i);
}

// Interpolate vector data inside or near the cell
template<int dim>
Vec3 LinearInterpolator<dim>::interp_vector(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

//    SimpleCell<dim> scell = cells[cell];
    array<unsigned, dim> scell = solution_indices[cell];

    // calculate weights or barycentric coordinates
    array<double,dim> bcc;
    if (c >= 0) bcc = get_bcc(Vec3(point), cell);
    else bcc = get_weights(point, cell);

    // Interpolate electric field
    Vec3 vector_i(0.0);
    for (int i = 0; i < dim; ++i)
        vector_i += solutions[scell[i]].vector * bcc[i];

    return vector_i;
}

// Interpolate scalar data inside or near the cell
template<int dim>
double LinearInterpolator<dim>::interp_scalar(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

//    SimpleCell<dim> scell = cells[cell];
    array<unsigned, dim> scell = solution_indices[cell];

    // calculate weights or barycentric coordinates
    array<double,dim> bcc;
    if (c >= 0) bcc = get_bcc(Vec3(point), cell);
    else bcc = get_weights(point, cell);

    // Interpolate potential
    double scalar_i(0.0);
    for (int i = 0; i < dim; ++i)
        scalar_i += solutions[scell[i]].scalar * bcc[i];

    return scalar_i;
}

// Find the cell which contains the point or is the closest to it
template<int dim>
int LinearInterpolator<dim>::locate_cell(const Point3 &point, const int cell_guess) const {
    // Check the guessed cell
    Vec3 vec_point(point);
    if (point_in_cell(vec_point, cell_guess)) return cell_guess;

    const int n_cells = cells.size();
    const int n_nbor_layers = 6;  // amount of nearest neighbouring layers that are checked before the full search

    vector<vector<int>> nbors(n_nbor_layers);
    vector<bool> cell_checked(n_cells);

    cell_checked[cell_guess] = true;

    // Check all cells on the given neighbouring layer
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of neighbour list
        if (layer == 0)
            nbors[0] = neighbours[cell_guess];
        else {
            for (int nbor : nbors[layer-1])
                if (nbor >= 0)
                    nbors[layer].insert(nbors[layer].end(), neighbours[nbor].begin(), neighbours[nbor].end());
        }

        // check whether some of the unchecked neighbouring cells surround the point
        for (int cell : nbors[layer])
            if (cell >= 0 && !cell_checked[cell]) {
                if (point_in_cell(vec_point, cell))
                    return cell;
                else
                    cell_checked[cell] = true;
            }
    }

    // If no success, loop through all the cells
    double min_distance2 = 1e100;
    int min_index = 0;

    for (int cell = 0; cell < n_cells; ++cell) {
        // If correct cell is found, we're done
        if (!cell_checked[cell] && point_in_cell(vec_point, cell))
            return cell;

        // Otherwise look for the cell whose centroid is closest to the point
        else {
            const double distance2 = point.distance2(centroids[cell]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = cell;
            }
        }
    }

    // If no perfect cell found, return the best.
    // Indicate the imperfectness with the minus sign
    return -min_index;
}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
template<int dim>
bool LinearInterpolator<dim>::average_sharp_nodes(const vector<vector<unsigned>>& nborlist) {
    // loop through the tetrahedral nodes
    for (int i = 0; i < nborlist.size(); ++i) {
        if (nborlist[i].size() == 0) continue;

        Point3 tetnode = (*nodes)[i];
        Vec3 vec(0);
        double w_sum = 0;

        // tetnode new solution will be the weighed average of the solutions on neighbouring nodes
        for (unsigned v : nborlist[i]) {
            double w = exp(decay_factor * tetnode.distance((*nodes)[v]));
            w_sum += w;
            vec += solutions[v].vector * w;
        }

        if (w_sum > 0) {
            vec *= (1.0 / w_sum);
            solutions[i].vector = vec;
            solutions[i].norm = vec.norm();
        }
    }

    return false;
}

/* Return the mapping between tetrahedral & hexahedral mesh nodes,
 nodes & hexahedral elements and nodes & element's vertices  */
template<int dim>
void LinearInterpolator<dim>::get_maps(vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs,
        dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh, const double eps) {

    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    const int n_femocs_nodes = nodes->size();
    const int n_dealii_nodes = tria->n_used_vertices();
    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;

    vector<int> node2hex(n_dealii_nodes), node2vert(n_dealii_nodes);
    tet2hex = vector<int>(n_femocs_nodes, -1);

    typename dealii::Triangulation<3>::active_vertex_iterator vertex = tria->begin_active_vertex();
    // Loop through tetrahedral mesh vertices
    for (int i = 0; i < n_femocs_nodes && vertex != tria->end_vertex(); ++i)
        if ( (*nodes)[i].distance2(vertex->vertex()) < eps ) {
            tet2hex[i] = vertex->vertex_index();
            vertex++;
        }

    // Loop through the hexahedral mesh elements
    typename dealii::DoFHandler<3>::active_cell_iterator cell;
    for (cell = dofh->begin_active(); cell != dofh->end(); ++cell)
        // Loop through all the vertices in the element
        for (int i = 0; i < n_verts_per_elem; ++i) {
            node2hex[cell->vertex_index(i)] = cell->active_cell_index();
            node2vert[cell->vertex_index(i)] = i;
        }

    // Generate lists of hexahedra and hexahedra nodes where the tetrahedra nodes are located
    cell_indxs.reserve(n_femocs_nodes);
    vert_indxs.reserve(n_femocs_nodes);
    for (int n : tet2hex)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }
}

// Pick the suitable write function based on the file type
template<int dim>
void LinearInterpolator<dim>::write(const string &file_name) const {
    if (!MODES.WRITEFILE) return;

    const int n_nodes = nodes->stat.n_tetnode;
    expect(n_nodes, "Zero nodes detected!");
    string ftype = get_file_type(file_name);

    ofstream outfile;
    outfile.setf(std::ios::scientific);
    outfile.precision(6);

    if (ftype == "movie") outfile.open(file_name, ios_base::app);
    else outfile.open(file_name);
    require(outfile.is_open(), "Can't open a file " + file_name);

    if (ftype == "xyz" || ftype == "movie")
        write_xyz(outfile, n_nodes);
    else if (ftype == "vtk")
        write_vtk(outfile, n_nodes);
    else
        require(false, "Unsupported file type: " + ftype);

    outfile.close();
}

// Output node data in .xyz format
template<int dim>
void LinearInterpolator<dim>::write_xyz(ofstream& out, const int n_nodes) const {
    out << n_nodes << endl;
    out << "LinearInterpolator properties=id:I:1:pos:R:3:marker:I:1:force:R:3:" <<
            norm_label << ":R:1:" << scalar_label << ":R:1" << endl;

    for (int i = 0; i < n_nodes; ++i)
        out << i << " " << (*nodes)[i] << " " << nodes->get_marker(i) << " " << solutions[i] << endl;
}

// Output interpolation cell data in .vtk format
template<int dim>
void LinearInterpolator<dim>::write_vtk(ofstream& out, const int n_nodes) const {
    const int celltype = get_cell_type();
    const int n_cells = cells.size();

    out << "# vtk DataFile Version 3.0\n";
    out << "# Medium data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << (*nodes)[i] << "\n";

    // Output the vertex indices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << cells[i] << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";

    out << "\nPOINT_DATA " << n_nodes << "\n";

    // write node IDs
    out << "SCALARS node-ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << i << "\n";

    // write node markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << nodes->get_marker(i) << "\n";

    // write vector norm data
    out << "SCALARS " + norm_label + " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << solutions[i].norm << "\n";

    // write scalar data
    out << "SCALARS " + scalar_label + " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << solutions[i].scalar << "\n";

    out << "\nCELL_DATA " << n_cells << "\n";

    // write cell IDs
    out << "SCALARS cell-ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_cells; ++i)
        out << i << "\n";
}

/* ==================================================================
 *  =================== TetrahedronInterpolator ====================
 * ================================================================== */

// Initialize data vectors
TetrahedronInterpolator::TetrahedronInterpolator(const TetgenMesh* m) :
        LinearInterpolator<4>(m), elems(&m->elems) {}

// Print statistics about solution on node points
void TetrahedronInterpolator::print_statistics() const {
    if (!MODES.VERBOSE) return;

    const int n_atoms = solutions.size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        if (solutions[i].norm == 0) continue;
        double s = solutions[i].scalar;
        Vec3 v = solutions[i].vector;

        vec += v; rms_vec += v * v;
        scalar += s; rms_scalar += s * s;
        n_points++;
    }

    vec *= (1.0 / n_points);
    rms_vec = Vec3(sqrt(rms_vec.x), sqrt(rms_vec.y), sqrt(rms_vec.z)) * (1.0 / n_points);
    scalar = scalar / n_points;
    rms_scalar = sqrt(rms_scalar) / n_points;

    stringstream stream;
    stream << "mean vector: \t" << vec;
    stream << "\n   rms vector: \t" << rms_vec;
    stream << "\n  mean & rms scalar: " << scalar << "\t" << rms_scalar;

    write_verbose_msg(stream.str());
}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
bool TetrahedronInterpolator::average_sharp_nodes(const bool vacuum) {
    vector<vector<unsigned int>> vorocells;
    mesh->calc_pseudo_3D_vorocells(vorocells, vacuum);
    return LinearInterpolator<4>::average_sharp_nodes(vorocells);
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::Laplace<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = nodes->size();
    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");
    const double eps = 1e-5 * elems->stat.edgemin;

    // Precompute tetrahedra to make interpolation faster
    reserve(n_nodes);
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tetNode2hexNode, cell_indxs, vert_indxs;
    get_maps(tetNode2hexNode, cell_indxs, vert_indxs,
            fem->get_triangulation(), fem->get_dof_handler(), eps);

    vector<dealii::Tensor<1, 3>> fields = fem->get_efield(cell_indxs, vert_indxs); // get list of electric fields
    vector<double> pot = fem->get_potential(cell_indxs, vert_indxs); // get list of electric potentials

    require(fields.size() == pot.size(),
            "Mismatch of vector sizes: " + to_string(fields.size()) + ", " + to_string(pot.size()));

    unsigned i = 0;
    for (int n : tetNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < fields.size(), "Invalid index: " + to_string(i));
            dealii::Tensor<1, 3> field = fields[i]; // that step needed to avoid complaints from Valgrind
            append_solution( Solution(Vec3(field[0], field[1], field[2]), pot[i++]) );
        }

        // In case of non-common node, store solution with error value
        else
            append_solution(Solution(0));
    }

    // remove the spikes in the solution
    if (average_sharp_nodes(true))
        return true;

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::CurrentsAndHeatingStationary<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = nodes->size();
    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");
    const double eps = 1e-5 * elems->stat.edgemin;

    // Precompute tetrahedra to make interpolation faster
    reserve(n_nodes);
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tetNode2hexNode, cell_indxs, vert_indxs;
    get_maps(tetNode2hexNode, cell_indxs, vert_indxs,
            fem->get_triangulation(), fem->get_dof_handler(), eps);

    vector<dealii::Tensor<1, 3>> rho = fem->get_current(cell_indxs, vert_indxs); // get list of current densities
    vector<double> temperature = fem->get_temperature(cell_indxs, vert_indxs); // get list of temperatures

    require(rho.size() == temperature.size(),
            "Mismatch of vector sizes: " + to_string(rho.size()) + ", "
                    + to_string(temperature.size()));

    unsigned i = 0;
    for (int n : tetNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < rho.size(), "Invalid index: " + to_string(i));
            Vec3 current(rho[i][0], rho[i][1], rho[i][2]);
            append_solution(Solution(current, temperature[i++]));
        }

        // In case of non-common node, store solution with error value
        else
            append_solution(Solution(0));
    }

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::CurrentsAndHeating<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = nodes->size();
    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");
    const double eps = 1e-5 * elems->stat.edgemin;

    // Precompute tetrahedra to make interpolation faster
    reserve(n_nodes);
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tetNode2hexNode, cell_indxs, vert_indxs;
    get_maps(tetNode2hexNode, cell_indxs, vert_indxs,
            fem->get_triangulation(), fem->get_dof_handler_current(), eps);

    vector<dealii::Tensor<1, 3>> rho = fem->get_current(cell_indxs, vert_indxs); // get list of current densities
    vector<double> temperature = fem->get_temperature(cell_indxs, vert_indxs); // get list of temperatures

    require(rho.size() == temperature.size(),
            "Mismatch of vector sizes: " + to_string(rho.size()) + ", "
                    + to_string(temperature.size()));

    unsigned i = 0;
    for (int n : tetNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < rho.size(), "Invalid index: " + to_string(i));
            Vec3 current(rho[i][0], rho[i][1], rho[i][2]);
            append_solution(Solution(current, temperature[i++]));
        }

        // In case of non-common node, store solution with error value
        else
            append_solution(Solution(0));
    }

    return false;
}

// Reserve memory for pre-compute data
void TetrahedronInterpolator::reserve_precompute(const int N) {
    LinearInterpolator<4>::reserve_precompute(N);

    det0.clear();
    det1.clear();
    det2.clear();
    det3.clear();
    det4.clear();
    tet_not_valid.clear();

    det0.reserve(N);
    det1.reserve(N);
    det2.reserve(N);
    det3.reserve(N);
    det4.reserve(N);
    tet_not_valid.reserve(N);
}

int TetrahedronInterpolator::opposite_node(const int tet, const int vert) const {
    const int dim = 4;
    vector<SimpleHex> hexs;
    for (int i = 0; i < dim; ++i)
        if (i != vert) hexs.push_back(mesh->hexahedra[dim * tet + i]);

    for (int node0 : hexs[0]) {
        if (nodes->get_marker(node0) != TYPES.FACECENTROID) continue;
        for (int node1 : hexs[1]) {
            if (node1 != node0) continue;
            for (int node2 : hexs[2])
                if (node2 == node0)
                    return node0;
        }
    }

    return -1;
}

//int TetrahedronInterpolator::opposite_node_vol2(const int tet, const int node) {
//    vector<SimpleHex> other_hexs;
//    for (int i = 0; i < 4; ++i)
//        if (i != node) other_hexs.push_back(mesh->hexahedra[4 * tet + i]);
//
//    for (int other_node : other_hexs[0]) {
//        bool found = true;
//        for (int i = 1; i <= 2; ++i) {
//            SimpleHex hex2 = other_hexs[i];
//            if ((hex2 != other_node) || (hex == other_node)) {
//                found = false;
//                break;
//            }
//        }
//        if (found)
//            return other_node;
//    }
//    return -1;
//}

// Precompute the data to tetrahedra to make later bcc calculations faster
void TetrahedronInterpolator::precompute() {
    const int n_elems = elems->size();
    const int n_nodes = nodes->size();
    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    reserve_precompute(n_elems);
    double d0, d1, d2, d3, d4;

    // Store the constant for smoothing
    decay_factor = -1.0 / elems->stat.edgemax;

    // Store the coordinates of tetrahedra vertices
    for (int i = 0; i < n_nodes; ++i)
        vertices.push_back((*nodes)[i]);

    // Calculate tetrahedra neighbours
    for (int i = 0; i < n_elems; ++i)
        neighbours.push_back(elems->get_neighbours(i));

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroids.push_back(elems->get_centroid(i));

    // Store tetrahedra to become free to interpolate independently from mesh state
    for (int i = 0; i < n_elems; ++i)
        cells.push_back((*elems)[i]);

    // Store the indices of face centroids where the solution is located
    for (int i = 0; i < n_elems; ++i) {
        array<unsigned, 4> opposite_nodes;
        for (int node = 0; node < 4; ++node)
            opposite_nodes[node] = opposite_node(i, node);
        solution_indices.push_back(opposite_nodes);
    }

    /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
     * barycentric coordinate of tetrahedra using the relations below */
    for (SimpleElement se : *elems) {
        Vec3 v1 = nodes->get_vec(se[0]);
        Vec3 v2 = nodes->get_vec(se[1]);
        Vec3 v3 = nodes->get_vec(se[2]);
        Vec3 v4 = nodes->get_vec(se[3]);

        /* =====================================================================================
         * det0 = |x1 y1 z1 1|
                  |x2 y2 z2 1|
                  |x3 y3 z3 1|
                  |x4 y4 z4 1|  */
        d0 = determinant(v1, v2, v3, v4);
        tet_not_valid.push_back(fabs(d0) < epsilon);
        det0.push_back(1.0 / d0);

        /* =====================================================================================
         * det1 = |x  y  z  1| = + x * |y2 z2 1| - y * |x2 z2 1| + z * |x2 y2 1| - |x2 y2 z2|
                  |x2 y2 z2 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));
        d2 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.z, v3.z, v4.z));
        d3 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y));
        d4 = determinant(Vec3(v2.x, v3.x, v4.x), Vec3(v2.y, v3.y, v4.y), Vec3(v2.z, v3.z, v4.z));

        det1.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det2 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x  y  z  1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x3 y3 z3 1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));
        d2 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.z, v3.z, v4.z));
        d3 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y));
        d4 = determinant(Vec3(v1.x, v3.x, v4.x), Vec3(v1.y, v3.y, v4.y), Vec3(v1.z, v3.z, v4.z));

        det2.push_back(Vec4(-d1, d2, -d3, d4));

        /* =====================================================================================
         * det3 = |x1 y1 z1 1| = + x * |y1 z1 1| - y * |x1 z1 1| + z * |x1 y1 1| - |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x  y  z  1|         |y4 z4 1|       |x4 z4 1|       |x4 y4 1|   |x4 y4 z4|
                  |x4 y4 z4 1|  */
        d1 = determinant(Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));
        d2 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.z, v2.z, v4.z));
        d3 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y));
        d4 = determinant(Vec3(v1.x, v2.x, v4.x), Vec3(v1.y, v2.y, v4.y), Vec3(v1.z, v2.z, v4.z));

        det3.push_back(Vec4(d1, -d2, d3, -d4));

        /* =====================================================================================
         * det4 = |x1 y1 z1 1| = - x * |y1 z1 1| + y * |x1 z1 1| - z * |x1 y1 1| + |x1 y1 z1|
                  |x2 y2 z2 1|         |y2 z2 1|       |x2 z2 1|       |x2 y2 1|   |x2 y2 z2|
                  |x3 y3 z3 1|         |y3 z3 1|       |x3 z3 1|       |x3 y3 1|   |x3 y3 z3|
                  |x  y  z  1|  */
        d1 = determinant(Vec3(v1.y, v2.y, v3.y), Vec3(v1.z, v2.z, v3.z));
        d2 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.z, v2.z, v3.z));
        d3 = determinant(Vec3(v1.x, v2.x, v3.x), Vec3(v1.y, v2.y, v3.y));
        d4 = determinant(v1, v2, v3);

        det4.push_back(Vec4(-d1, d2, -d3, d4));
    }
}

// Check with barycentric coordinates whether the point is inside the i-th tetrahedron
bool TetrahedronInterpolator::point_in_cell(const Vec3 &point, const int i) const {
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

    // Ignore co-planar tetrahedra
    // no need to check because Tetgen guarantees non-co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    const Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < 0) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < 0) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < 0) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < 0) return false;

    // All bcc-s are >= 0, so point is inside the tetrahedron
    return true;
}

// Calculate barycentric coordinates for point with respect to i-th tetrahedron
array<double,4> TetrahedronInterpolator::get_bcc(const Vec3& point, const int i) const {
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

    const Vec4 pt(point, 1);
    const double bcc1 = det0[i] * pt.dotProduct(det1[i]);
    const double bcc2 = det0[i] * pt.dotProduct(det2[i]);
    const double bcc3 = det0[i] * pt.dotProduct(det3[i]);
    const double bcc4 = det0[i] * pt.dotProduct(det4[i]);

//    return array<double,4> {bcc1, bcc2, bcc3, bcc4};
    return array<double,4> {1.0 - 3.0 * bcc1, 1.0 - 3.0 * bcc2, 1.0 - 3.0 * bcc3, 1.0 - 3.0 * bcc4};
}

/* ==================================================================
 *  ===================== TriangleInterpolator ======================
 * ================================================================== */

// Initialize data vectors
TriangleInterpolator::TriangleInterpolator(const TetgenMesh* m) :
        LinearInterpolator<3>(m), faces(&m->faces) {}

// leave only the solution in the nodes and centroids of triangles
bool TriangleInterpolator::clean_nodes() {
    const int n_nodes = nodes->size();

    vector<bool> node_not_in_quads(n_nodes, true);
    for (SimpleQuad quad : mesh->quads)
        for (int node : quad)
            node_not_in_quads[node] = false;

    for (int node = 0; node < n_nodes; ++node)
        if (node_not_in_quads[node])
            solutions[node] = Solution(0);

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TriangleInterpolator::extract_solution(fch::Laplace<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = nodes->size();
    expect(n_nodes > 0, "Interpolator expects non-empty mesh!");
    const double eps = 1e-5 * faces->stat.edgemin;

    // Precompute triangles to make interpolation faster
    reserve(n_nodes);
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> triNode2hexNode, cell_indxs, vert_indxs;
    get_maps(triNode2hexNode, cell_indxs, vert_indxs, fem->get_triangulation(), fem->get_dof_handler(), eps);

    vector<dealii::Tensor<1, 3>> fields = fem->get_efield(cell_indxs, vert_indxs); // get list of electric fields
    vector<double> pot = fem->get_potential(cell_indxs, vert_indxs); // get list of electric potentials

    require(fields.size() == pot.size(),
            "Mismatch of vector sizes: " + to_string(fields.size()) + ", " + to_string(pot.size()));

    unsigned i = 0;
    for (int n : triNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < fields.size(), "Invalid index: " + to_string(i));
            dealii::Tensor<1, 3> field = fields[i]; // that step needed to avoid complaints from Valgrind
            append_solution( Solution(Vec3(field[0], field[1], field[2]), pot[i++]) );
        }

        // In case of non-common node, store solution with error value
        else
            append_solution(Solution(0));
    }

    // remove the spikes in the solution
    if (average_sharp_nodes(true))
        return true;

    // leave only the solution in the vertices and centroids of triangles
    return clean_nodes();
}

// Interpolate scalar value inside or near cell
//  whose total sum must remain conserved after interpolations
void TriangleInterpolator::interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) {
    const int n_atoms = atoms.size();
    const int n_nodes = nodes->size();

    vector<double> bcc_sum(n_nodes); // sum of barycentric coordinates from given node
    vector<int> atom2cell(n_atoms);  // map storing the face indices that correspond to atom sequence
    scalars = vector<double>(n_atoms);
    array<double,3> bcc;

    // calculate the sum of all the weights in all the nodes
    // it is neccesary to ensure that the total sum of a interpolated scalar does not change
    int cell = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms[i].point;
        // Find the cell that matches best to the point
        cell = abs(locate_cell(point, cell));
        // calculate barycentric coordinates
        bcc = get_bcc(Vec3(point), cell);
        // store the cell to make next step faster
        atom2cell[i] = cell;
        // append barycentric weight to the sum of all weights from given node
        int j = 0;
        for (int node : (*faces)[cell])
            bcc_sum[node] += bcc[j++];
    }

    // force bcc_sum in the location of unused nodes to some non-zero value
    // to avoid nan-s in bcc[i]/bcc_sum[i]
    for (int i = 0; i < n_nodes; ++i)
        if (bcc_sum[i] == 0)
            bcc_sum[i] = 1;

    // perform actual interpolation for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        int cell = atom2cell[i];
        // calculate barycentric coordinates
        bcc = get_bcc(Vec3(atoms[i].point), cell);
        // perform interpolation
        int j = 0;
        for (int node : (*faces)[cell])
            scalars[i] += solutions[node].scalar * bcc[j++] / bcc_sum[node];
    }
}

// Calculate charges on surface faces using direct solution in the face centroids
void TriangleInterpolator::calc_charges(const double E0) {
    const int n_quads_per_triangle = 3;
    const int n_faces = faces->size();
    const int n_nodes = nodes->size();
    const double sign = fabs(E0) / E0;

    vector<double> charges(n_nodes), areas(n_nodes);

    // create triangle index to its centroid index mapping
    vector<int> tri2centroid(n_faces);
    for (int face = 0; face < n_faces; ++face) {
        for (int node : mesh->quads[n_quads_per_triangle * face])
            if (nodes->get_marker(node) == TYPES.FACECENTROID) {
                tri2centroid[face] = node;
                break;
            }
    }

    // find the nodes on the surface perimeter that are each-other's periodic images
    vector<int> periodic_nodes(n_nodes, -1);
    vector<bool> node_on_edge(n_nodes);
    for (SimpleEdge edge : mesh->edges)
        for (int node : edge)
            node_on_edge[node] = true;

    for (int i = 0; i < vector_sum(node_on_edge); i += 2) {
        periodic_nodes[i] = i+1;
        periodic_nodes[i+1] = i;
    }

    // Calculate the charges and areas in the centroids and vertices of triangles
    for (int face = 0; face < n_faces; ++face) {
        int centroid_indx = tri2centroid[face];
        double area = faces->get_area(face);
        double charge = eps0 * sign * area * solutions[centroid_indx].norm;
        solutions[centroid_indx].norm = area;
        solutions[centroid_indx].scalar = charge;

        charge *= 1.0 / n_quads_per_triangle;
        area *= 1.0 / n_quads_per_triangle;

        // the charge on triangular node is the sum of the charges of quadrangules that surround the node
        for (int node : (*faces)[face]) {
            const int periodic_node = periodic_nodes[node];
            charges[node] += charge;
            areas[node] += area;
            if (periodic_node >= 0) {
                charges[periodic_node] += charge;
                areas[periodic_node] += area;
            }
        }
    }

    // transfer charges and areas to solutions vector
    for (SimpleFace face : *faces) {
        for (int node : face) {
            solutions[node].norm = areas[node];
            solutions[node].scalar = charges[node];
        }
    }
}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
bool TriangleInterpolator::average_sharp_nodes(const bool vacuum) {
    vector<vector<unsigned int>> vorocells;
    mesh->calc_pseudo_3D_vorocells(vorocells, vacuum);
    return LinearInterpolator<3>::average_sharp_nodes(vorocells);
}

// Reserve memory for precompute data
void TriangleInterpolator::reserve_precompute(const int n) {
    LinearInterpolator<3>::reserve_precompute(n);

    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear();  pvec.reserve(n);
    norms.clear(); norms.reserve(n);
    max_distance.clear(); max_distance.reserve(n);
}

// Precompute the data needed to calculate the distance of points from surface in the direction of triangle norms
void TriangleInterpolator::precompute() {
    const int n_faces = faces->size();
    const int n_nodes = nodes->size();

    // Reserve memory for precomputation data
    reserve_precompute(n_faces);

    // Store the constant for smoothing
    decay_factor = -1.0 / faces->stat.edgemax;

    // Store the coordinates of tetrahedral vertices
    // Tetrahedral not triangular, because the only solid knowledge about the triangle nodes
    // is that they are among tetrahedral nodes
    for (int i = 0; i < n_nodes; ++i)
        vertices.push_back((*nodes)[i]);

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace sface = (*faces)[i];

        Vec3 v0 = nodes->get_vec(sface[0]);
        Vec3 v1 = nodes->get_vec(sface[1]);
        Vec3 v2 = nodes->get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = faces->get_norm(i).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);

        // store faces and norms to be able to interpolate independently of mesh state
        cells.push_back((*faces)[i]);
        norms.push_back(faces->get_norm(i));
        max_distance.push_back(e2.norm());
        centroids.push_back(faces->get_centroid(i));

        // calculate the neighbour list for triangles
        for (int j = i+1; j < n_faces; ++j)
            if (sface.edge_neighbor((*faces)[j])) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }

        // Store the indices of edge centroids where the solution is located
        array<unsigned, 3> opposite_nodes;
        for (int node = 0; node < 3; ++node)
            opposite_nodes[node] = opposite_node(i, node);
        solution_indices.push_back(opposite_nodes);
    }
}

int TriangleInterpolator::near_surface(const Vec3& point, const double r_cut) const {
    require(r_cut > 0, "Invalid distance from surface: " + to_string(r_cut));

    for (int face = 0; face < cells.size(); ++face) {
        const double dist = distance(point, face);
        if (dist >= -0.3*r_cut && dist <= r_cut) return face;
    }

    return -1;
}

int TriangleInterpolator::opposite_node(const int tri, const int vert) const {
    return (*faces)[tri][vert];
//    const int dim = 3;
//    vector<SimpleQuad> quads;
//    for (int i = 0; i < dim; ++i)
//        if (i != vert) quads.push_back(mesh->quads[dim * tri + i]);
//
//    for (int node0 : quads[0]) {
//        if (nodes->get_marker(node0) != TYPES.EDGECENTROID) continue;
//        for (int node1 : quads[1])
//            if (node1 == node0)
//                return node0;
//    }
//
//    return -1;
}

Vec3 TriangleInterpolator::get_norm(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return norms[i];
}

double TriangleInterpolator::fast_distance(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    return edge2[face].dotProduct(qvec);
}

double TriangleInterpolator::distance(const Vec3& point, const int face) const {
    // Constants to specify the tolerances in searching outside the triangle
    const static double zero = -0.1;
    const static double one = 1.1;

    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < zero || u > one) return 1e100;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = norms[face].dotProduct(qvec);
    if (v < zero || u + v > one) return 1e100; // Check second & third barycentric coordinate

    // return the distance from point to triangle
    return edge2[face].dotProduct(qvec);
}

void TriangleInterpolator::write_vtk(ofstream& out, const int n_nodes) const {
    LinearInterpolator<3>::write_vtk(out, n_nodes);

    // write face norms
    out << "VECTORS norm double\n";
    for (int i = 0; i < cells.size(); ++i)
        out << norms[i] << "\n";
}

// Calculate barycentric coordinates for a projection of a point inside the triangle
array<double,3> TriangleInterpolator::get_bcc(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double u = tvec.dotProduct(pvec[face]);
    double v = qvec.dotProduct(norms[face]);
    double w = 1.0 - u - v;

    return array<double,3> {w, u, v};
//    return array<double,3> {1.0 - 2.0 * w, 1.0 - 2.0 * u, 1.0 - 2.0 * v};
}

// Check whether the projection of a point is inside the triangle
// It is separate routine from get_bcc to achieve better performance
bool TriangleInterpolator::point_in_cell(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < 0 || u > 1) return false;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = qvec.dotProduct(norms[face]);
    if (v < 0 || u + v > 1) return false; // Check second & third barycentric coordinate

    // finally check the distance of the point from the triangle
    return fabs(qvec.dotProduct(edge2[face])) < max_distance[face];
}

void TriangleInterpolator::print_statistics(const double Q) {
    if (!MODES.VERBOSE) return;
    if (nodes->size() == size()) {
        expect(false, "Mismatch between interpolator and mesh sizes: " + to_string(size()) + ", "
                + to_string(nodes->size()));
        return;
    }

    double q_sum = 0;
    for (int i = 0; i < size(); ++i) {
        const double scalar = get_scalar(i);
        if (nodes->get_marker(i) == TYPES.TETNODE)
            q_sum += scalar;
    }

    stringstream stream; stream << fixed << setprecision(3);
    stream << "Q / sum(" << scalar_label << ") = " << Q << " / " << q_sum << " = " << Q/q_sum;
    write_verbose_msg(stream.str());
}

template class LinearInterpolator<3> ;
template class LinearInterpolator<4> ;

} // namespace femocs
