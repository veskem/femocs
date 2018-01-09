/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "../include/Interpolator.h"

#include "Macros.h"
#include "float.h"

using namespace std;
namespace femocs {

/* ==================================================================
 *  ======================== Interpolator ==========================
 * ================================================================== */

// Pick the suitable write function based on the file type
void Interpolator::write(const string &file_name) const {
    if (!MODES.WRITEFILE) return;

    ofstream outfile;
    outfile.setf(std::ios::scientific);
    outfile.precision(6);

    string ftype = get_file_type(file_name);
    if (ftype == "movie") outfile.open(file_name, ios_base::app);
    else outfile.open(file_name);
    require(outfile.is_open(), "Can't open a file " + file_name);

    if (ftype == "xyz" || ftype == "movie")
        write_xyz(outfile);
    else if (ftype == "vtk")
        write_vtk(outfile);
    else
        require(false, "Unsupported file type: " + ftype);

    outfile.close();
}

void Interpolator::write_xyz(ofstream& out) const {
    const int n_nodes = nodes->size(); // nodes->stat.n_tetnode;
    expect(n_nodes, "Zero nodes detected!");

    out << n_nodes << endl;
    out << "Interpolator properties=id:I:1:pos:R:3:marker:I:1:force:R:3:" <<
            norm_label << ":R:1:" << scalar_label << ":R:1" << endl;

    for (int i = 0; i < n_nodes; ++i)
        out << i << " " << (*nodes)[i] << " " << nodes->get_marker(i) << " " << solutions[i] << endl;
}

void Interpolator::write_vtk(ofstream& out) const {
    const int n_nodes = nodes->size(); // nodes->stat.n_tetnode;
    expect(n_nodes, "Zero nodes detected!");

    const int celltype = get_cell_type();
    const int n_cells = neighbours.size();

    out << "# vtk DataFile Version 3.0\n";
    out << "# Medium data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << (*nodes)[i] << "\n";

    // Output the vertex indices and cell types
    write_cells(out);

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

int Interpolator::locate_cell(const Point3 &point, const int cell_guess) const {
    vector<bool> cell_checked(neighbours.size());
    return locate_cell(point, cell_guess, cell_checked);
}

// Find the cell which contains the point or is the closest to it
int Interpolator::locate_cell(const Point3 &point, const int cell_guess, vector<bool>& cell_checked) const {
    // Check the guessed cell
    Vec3 vec_point(point);
    if (point_in_cell(vec_point, cell_guess)) return cell_guess;
    cell_checked[cell_guess] = true;

    const int n_cells = neighbours.size();
    const int n_nbor_layers = 6;  // amount of nearest neighbouring layers that are checked before the full search
    vector<vector<int>> nbors(n_nbor_layers);

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
bool Interpolator::average_sharp_nodes(const vector<vector<unsigned>>& nborlist) {
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

/* Calculate mapping between Femocs & deal.II mesh nodes,
 nodes & hexahedral elements and nodes & element's vertices */
void Interpolator::get_maps(vector<int>& femocs2deal, vector<int>& cell_indxs, vector<int>& vert_indxs,
        dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh) const {

    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;

    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    const int n_femocs_nodes = nodes->size();
    expect(n_femocs_nodes > 0, "Interpolator expects non-empty mesh!");
    const int n_dealii_nodes = tria->n_used_vertices();

    vector<int> node2hex(n_dealii_nodes), node2vert(n_dealii_nodes);
    femocs2deal = vector<int>(n_femocs_nodes, -1);

    typename dealii::Triangulation<3>::active_vertex_iterator vertex = tria->begin_active_vertex();
    // Loop through tetrahedral mesh vertices
    for (int i = 0; i < n_femocs_nodes && vertex != tria->end_vertex(); ++i)
        if ( (*nodes)[i].distance2(vertex->vertex()) < vertex_epsilon ) {
            femocs2deal[i] = vertex->vertex_index();
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
    for (int n : femocs2deal)
        if (n >= 0) {
            cell_indxs.push_back(node2hex[n]);
            vert_indxs.push_back(node2vert[n]);
        }
}

void Interpolator::store_solution(const vector<int>& femocs2deal,
        const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data) {

    require(vec_data.size() == scal_data.size(), "Mismatch of vector sizes: "
            + to_string(vec_data.size()) + ", " + to_string(scal_data.size()));

    // Reserve memory for data
    reserve(nodes->size());

    unsigned i = 0;
    for (int n : femocs2deal) {
        // If there is a common node between Femocs and deal.II meshes, store actual solution
        if (n >= 0) {
            require(i < vec_data.size(), "Invalid index: " + to_string(i));
            dealii::Tensor<1, 3> vec = vec_data[i]; // that step needed to avoid complaints from Valgrind
            append_solution( Solution(Vec3(vec[0], vec[1], vec[2]), scal_data[i++]) );
        }

        // In case of non-common node, store solution with error value
        else
            append_solution(Solution(0));
    }
}

bool Interpolator::extract_solution(fch::Laplace<3>* fem) {
    require(fem, "NULL pointer can't be handled!");

    // Precompute cells to make interpolation faster
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem->get_triangulation(), fem->get_dof_handler());

    // Read and store the electric field and potential from FEM solver
    store_solution(femocs2deal, fem->get_efield(cell_indxs, vert_indxs),
            fem->get_potential(cell_indxs, vert_indxs));

    // Remove the spikes from the solution
    return average_sharp_nodes(true);
}

bool Interpolator::extract_solution(fch::CurrentsAndHeatingStationary<3>* fem) {
    require(fem, "NULL pointer can't be handled!");

    // Precompute cells to make interpolation faster
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem->get_triangulation(), fem->get_dof_handler());

    // Read and store current densities and temperatures from FEM solver
    store_solution(femocs2deal, fem->get_current(cell_indxs, vert_indxs),
            fem->get_temperature(cell_indxs, vert_indxs));

    return false;
}

bool Interpolator::extract_solution(fch::CurrentsAndHeating<3>* fem) {
    require(fem, "NULL pointer can't be handled!");

    // Precompute cells to make interpolation faster
    precompute();

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> femocs2deal, cell_indxs, vert_indxs;
    get_maps(femocs2deal, cell_indxs, vert_indxs, fem->get_triangulation(), fem->get_dof_handler_current());

    // Read and store current densities and temperatures from FEM solver
    store_solution(femocs2deal, fem->get_current(cell_indxs, vert_indxs),
            fem->get_temperature(cell_indxs, vert_indxs));

    return false;
}

/* ==================================================================
 *  =================== TetrahedronInterpolator ====================
 * ================================================================== */

// Print statistics about solution on node points
template<int rank>
void TetrahedronInterpolator<rank>::print_statistics() const {
    if (!MODES.VERBOSE) return;

    const int n_atoms = Interpolator::size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        if ( Interpolator::get_vector_norm(i) == 0) continue;
        double s = Interpolator::get_scalar(i);
        Vec3 v = Interpolator::get_vector(i);

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
template<int rank>
bool TetrahedronInterpolator<rank>::average_sharp_nodes(const bool vacuum) {
    vector<vector<unsigned int>> vorocells;
    this->mesh->calc_pseudo_3D_vorocells(vorocells, vacuum);
    return Interpolator::average_sharp_nodes(vorocells);
}

// Reserve memory for pre-compute data
template<int rank>
void TetrahedronInterpolator<rank>::reserve_precompute(const int N) {
    Interpolator::reserve_precompute(N);
    cells.clear();
    cells.reserve(N);

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

// Precompute the data to tetrahedra to make later bcc calculations faster
template<int rank>
void TetrahedronInterpolator<rank>::precompute() {
    const int n_elems = elems->size();
    const int n_nodes = this->nodes->size();

    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    double d0, d1, d2, d3, d4;

    reserve_precompute(n_elems);

    // Store the constant for smoothing
    this->decay_factor = -1.0 / elems->stat.edgemax;

    // Store max allowed distance between identical vertices
    this->vertex_epsilon = 1e-5 * elems->stat.edgemin;

    // Store the coordinates of tetrahedra vertices
    for (int i = 0; i < n_nodes; ++i)
        this->vertices.push_back((*(this->nodes))[i]);

    for (int i = 0; i < n_elems; ++i) {
        SimpleElement se = (*elems)[i];

        // Calculate tetrahedra neighbours
        this->neighbours.push_back(this->elems->get_neighbours(i));
        // Calculate centroids of tetrahedra
        this->centroids.push_back(this->elems->get_centroid(i));
        // Store tetrahedra to get rid of mesh dependency
        this->cells.push_back(get_cell(i));

        /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
         * barycentric coordinate of tetrahedra using the relations below */

        Vec3 v1 = this->nodes->get_vec(se[0]);
        Vec3 v2 = this->nodes->get_vec(se[1]);
        Vec3 v3 = this->nodes->get_vec(se[2]);
        Vec3 v4 = this->nodes->get_vec(se[3]);

        /* =====================================================================================
         * det0 = |x1 y1 z1 1|
                  |x2 y2 z2 1|
                  |x3 y3 z3 1|
                  |x4 y4 z4 1|  */
        d0 = determinant(v1, v2, v3, v4);
        tet_not_valid.push_back(fabs(d0) < this->coplanar_epsilon);
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
template<int rank>
bool TetrahedronInterpolator<rank>::point_in_cell(const Vec3 &point, const int i) const {
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

// Calculate distance-dependent weights for a point with respect to the cell
template<int rank>
void TetrahedronInterpolator<rank>::get_weights(array<double,rank>& weights, const Point3 &point, const SimpleCell<rank>& scell) const {
    double w_sum = 0;
    for (int i = 0; i < rank; ++i) {
        double w = exp(decay_factor * point.distance(vertices[scell[i]]));
        weights[i] = w;
        w_sum += w;
    }

    // check for zero and nan
    require(w_sum > 0 && w_sum == w_sum, "Invalid interpolation weight: " + to_string(w_sum));
    w_sum = 1.0 / w_sum;
    for (int i = 0; i < rank; ++i)
        weights[i] *= w_sum;
}

// Interpolate both scalar and vector data inside or near the cell
template<int rank>
Solution TetrahedronInterpolator<rank>::interp_solution(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

    SimpleCell<rank> scell = cells[cell];

    // calculate weights or barycentric coordinates
    array<double,rank> weights;
    if (c >= 0) get_shape_functions(weights, Vec3(point), cell);
    else get_weights(weights, point, scell);

    // Interpolate electric field
    Vec3 vector_i(0.0);
    for (int i = 0; i < rank; ++i)
        vector_i += solutions[scell[i]].vector * weights[i];

    // Interpolate potential
    double scalar_i(0.0);
    for (int i = 0; i < rank; ++i)
        scalar_i += solutions[scell[i]].scalar * weights[i];

    return Solution(vector_i, scalar_i);
}

/* ==================================================================
 *  ===================== LinTetInterpolator =======================
 * ================================================================== */

void LinTetInterpolator::get_shape_functions(array<double,4>& sf, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < det0.size(), "Index out of bounds: " + to_string(tet));

    const Vec4 pt(point, 1);
    sf[0] = det0[tet] * pt.dotProduct(det1[tet]);
    sf[1] = det0[tet] * pt.dotProduct(det2[tet]);
    sf[2] = det0[tet] * pt.dotProduct(det3[tet]);
    sf[3] = det0[tet] * pt.dotProduct(det4[tet]);
}

SimpleCell<4> LinTetInterpolator::get_cell(const int i) const {
    if (elems->size() <= i) return SimpleElement(0);
    return (*elems)[i];
}

/* ==================================================================
 *  ===================== QuadTetInterpolator ======================
 * ================================================================== */

void QuadTetInterpolator::get_shape_functions(array<double,10>& sf, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < det0.size(), "Index out of bounds: " + to_string(tet));

    const Vec4 pt(point, 1);
    const double b1 = det0[tet] * pt.dotProduct(det1[tet]);
    const double b2 = det0[tet] * pt.dotProduct(det2[tet]);
    const double b3 = det0[tet] * pt.dotProduct(det3[tet]);
    const double b4 = det0[tet] * pt.dotProduct(det4[tet]);

    sf[0] =  b1 * (2 * b1 - 1);
    sf[1] =  b2 * (2 * b2 - 1);
    sf[2] =  b3 * (2 * b3 - 1);
    sf[3] =  b4 * (2 * b4 - 1);
    sf[4] =  4 * b1 * b2;
    sf[5] =  4 * b2 * b3;
    sf[6] =  4 * b3 * b1;
    sf[7] =  4 * b1 * b4;
    sf[8] =  4 * b2 * b4;
    sf[9] =  4 * b3 * b4;

//    sf = {b1, b2, b3, b4, 0, 0, 0, 0, 0, 0};
}

SimpleCell<10> QuadTetInterpolator::get_cell(const int tet) const {
    if (mesh->hexahedra.size() <= tet) 
        return QuadraticTet(0);
    
    const int n_hexs_per_tet = 4;
    array<vector<unsigned>,4> edge_nodes;

    // locate hexahedral nodes that are located in the middle of edges
    for (int i = 0; i < n_hexs_per_tet; ++i) {
        for (int hexnode : mesh->hexahedra[n_hexs_per_tet * tet + i])
            if (nodes->get_marker(hexnode) == TYPES.EDGECENTROID)
                edge_nodes[i].push_back(hexnode);
    }

    // find second order nodes
    const int n5 = common_entry(edge_nodes[0], edge_nodes[1]);
    const int n6 = common_entry(edge_nodes[1], edge_nodes[2]);
    const int n7 = common_entry(edge_nodes[2], edge_nodes[0]);
    const int n8 = common_entry(edge_nodes[0], edge_nodes[3]);
    const int n9 = common_entry(edge_nodes[1], edge_nodes[3]);
    const int n10= common_entry(edge_nodes[2], edge_nodes[3]);

    return QuadraticTet((*elems)[tet], n5, n6, n7, n8, n9, n10);
}

bool QuadTetInterpolator::average_and_check_sharp_nodes(const bool vacuum) {
    vector<vector<unsigned int>> vorocells;
    mesh->calc_pseudo_3D_vorocells(vorocells, vacuum);

    const int n_hexs = mesh->hexahedra.size();

    double total_before = 0;
    for (int i = 0; i < n_hexs; ++i)
        total_before += integrate(i);

    bool success = Interpolator::average_sharp_nodes(vorocells);

    double total_after = 0;
    int data_size = 0;
    for (int i = 0; i < n_hexs; ++i) {
        double energy = integrate(i);
        total_after += energy;
        if (energy != 0)
            data_size++;
    }

//    if (MODES.VERBOSE)
        printf("\n  E_before=%.3f, E_after=%.3f, diff=%.3f\%\n",
                total_before / data_size, total_after / data_size, 100.0 * (total_after / total_before - 1.0));

    return success;
}

void QuadTetInterpolator::hex_shape_function(array<double,8>& sf, const double u, const double v, const double w) const {
    const double n1 = (1-u) * (1-v) * (1-w) / 8;
    const double n2 = (1+u) * (1-v) * (1-w) / 8;
    const double n3 = (1+u) * (1+v) * (1-w) / 8;
    const double n4 = (1-u) * (1+v) * (1-w) / 8;
    const double n5 = (1-u) * (1-v) * (1+w) / 8;
    const double n6 = (1+u) * (1-v) * (1+w) / 8;
    const double n7 = (1+u) * (1+v) * (1+w) / 8;
    const double n8 = (1-u) * (1+v) * (1+w) / 8;
    sf = {n1, n2, n3, n4, n5, n6, n7, n8};
}

double QuadTetInterpolator::integrate(const int hex_index) const {
    const int n_nodes_per_hex = 8;
    const int n_steps = 10;
    const double step = 2.0 / n_steps;
    const SimpleHex elem = mesh->hexahedra[hex_index];

    array<double,8> sf;
    array<Vec3,8> nodal_field;
    for (int i = 0; i < n_nodes_per_hex; ++i)
        nodal_field[i] = solutions[elem[i]].vector;

    double total_energy = 0.0;

    for (double u = -1; u < 1; u += step)
        for (double v = -1; v < 1; v += step)
            for (double w = -1; w < 1; w += step) {
                hex_shape_function(sf, u, v, w);
                Vec3 dEnergy(0);
                for (int i = 0; i < n_nodes_per_hex; ++i)
                    dEnergy += nodal_field[i] * sf[i];

                total_energy += dEnergy.norm2();
            }

    return step * step * step * total_energy / 2.0 / eps0;
}

// Locate the index of node that is in the centroid of opposite face of the given tetrahedral vertex
int QuadTetInterpolator::opposite_node(const int tet, const int vert) const {
    const int n_hexs_per_tet = 4;
    vector<SimpleHex> hexs;
    for (int i = 0; i < n_hexs_per_tet; ++i)
        if (i != vert) hexs.push_back(mesh->hexahedra[n_hexs_per_tet * tet + i]);

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

/* ==================================================================
 *  ===================== TriangleInterpolator =====================
 * ================================================================== */

// Interpolate scalar value inside or near cell
//  whose total sum must remain conserved after interpolations
template<int rank>
void TriangleInterpolator<rank>::interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) const {
    const int n_atoms = atoms.size();
    const int n_nodes = this->nodes->size();

    vector<double> bcc_sum(n_nodes); // sum of barycentric coordinates from given node
    vector<int> atom2cell(n_atoms);  // map storing the face indices that correspond to atom sequence
    scalars = vector<double>(n_atoms);
    array<double,rank> weights;

    // calculate the sum of all the weights in all the nodes
    // it is neccesary to ensure that the total sum of a interpolated scalar does not change
    int cell = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms[i].point;
        // Find the cell that matches best to the point
        cell = abs(this->locate_cell(point, cell));
        // calculate barycentric coordinates
        get_shape_functions(weights, Vec3(point), cell);
        // store the cell to make next step faster
        atom2cell[i] = cell;
        // append barycentric weight to the sum of all weights from given node
        int j = 0;
        for (int node : get_cell(cell))
            bcc_sum[node] += weights[j++];
    }

    // force bcc_sum in the location of unused nodes to some non-zero value
    // to avoid nan-s in weights[i]/bcc_sum[i]
    for (int i = 0; i < n_nodes; ++i)
        if (bcc_sum[i] == 0)
            bcc_sum[i] = 1;

    // perform actual interpolation for all the atoms
    for (int i = 0; i < n_atoms; ++i) {
        int cell = atom2cell[i];
        // calculate barycentric coordinates
        get_shape_functions(weights, Vec3(atoms[i].point), cell);
        // perform interpolation
        int j = 0;
        for (int node : get_cell(cell))
            scalars[i] += this->get_scalar(node) * weights[j++] / bcc_sum[node];
    }
}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
template<int rank>
bool TriangleInterpolator<rank>::average_sharp_nodes(const bool vacuum) {
    vector<vector<unsigned int>> vorocells;
    this->mesh->calc_pseudo_3D_vorocells(vorocells, vacuum);
    return Interpolator::average_sharp_nodes(vorocells);
}

// Reserve memory for precompute data
template<int rank>
void TriangleInterpolator<rank>::reserve_precompute(const int n) {
    Interpolator::reserve_precompute(n);
    cells.clear();
    cells.reserve(n);

    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear();  pvec.reserve(n);
    norms.clear(); norms.reserve(n);
    max_distance.clear(); max_distance.reserve(n);
}

// Precompute the data needed to calculate the distance of points from surface in the direction of triangle norms
template<int rank>
void TriangleInterpolator<rank>::precompute() {
    const int n_faces = faces->size();
    const int n_nodes = this->nodes->size();

    // Reserve memory for precomputation data
    reserve_precompute(n_faces);

    // Store the constant for smoothing
    this->decay_factor = -1.0 / faces->stat.edgemax;

    // Store max allowed distance between identical vertices
    this->vertex_epsilon = 1e-5 * faces->stat.edgemin;

    // Store the coordinates of tetrahedral vertices
    // Tetrahedral not triangular, because the only solid knowledge about the triangle nodes
    // is that they are among tetrahedral nodes
    for (int i = 0; i < n_nodes; ++i)
        this->vertices.push_back( (*(this->nodes))[i] );

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace sface = (*faces)[i];

        Vec3 v0 = this->nodes->get_vec(sface[0]);
        Vec3 v1 = this->nodes->get_vec(sface[1]);
        Vec3 v2 = this->nodes->get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = faces->get_norm(i).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);

        // store triangles to get rid of dependence on mesh state
        this->cells.push_back(get_cell(i));
        // calculate norms of triangles
        this->norms.push_back(faces->get_norm(i));
        // store max distance from given triangle
        this->max_distance.push_back(e2.norm());
        // calculate centroids of triangles
        this->centroids.push_back(faces->get_centroid(i));

        // calculate the neighbour list for triangles
        for (int j = i+1; j < n_faces; ++j)
            if (sface.edge_neighbor((*faces)[j])) {
                this->neighbours[i].push_back(j);
                this->neighbours[j].push_back(i);
            }
    }
}

template<int rank>
int TriangleInterpolator<rank>::near_surface(const Vec3& point, const double r_cut) const {
    require(r_cut > 0, "Invalid distance from surface: " + to_string(r_cut));

    for (int face = 0; face < this->cells.size(); ++face) {
        const double dist = distance(point, face);
        if (dist >= -0.3*r_cut && dist <= r_cut) return face;
    }

    return -1;
}

template<int rank>
double TriangleInterpolator<rank>::fast_distance(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    return edge2[face].dotProduct(qvec);
}

template<int rank>
double TriangleInterpolator<rank>::distance(const Vec3& point, const int face) const {
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

template<int rank>
void TriangleInterpolator<rank>::write_vtk(ofstream& out) const {
    Interpolator::write_vtk(out);

    // write face norms
    out << "VECTORS norm double\n";
    for (int i = 0; i < this->cells.size(); ++i)
        out << norms[i] << "\n";
}

// Check whether the projection of a point is inside the triangle
// It is separate routine from get_bcc to achieve better performance
template<int rank>
bool TriangleInterpolator<rank>::point_in_cell(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < 0 || u > 1) return false;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = qvec.dotProduct(norms[face]);
    if (v < 0 || u + v > 1) return false; // Check second & third barycentric coordinate

    // finally check the distance of the point from the triangle
    return fabs(qvec.dotProduct(edge2[face])) < max_distance[face];
}

template<int rank>
void TriangleInterpolator<rank>::print_statistics(const double Q) {
    if (!MODES.VERBOSE) return;
    const int n_nodes = this->size();
    if (this->nodes->size() == n_nodes) {
        expect(false, "Mismatch between interpolator and mesh sizes: " + to_string(n_nodes) + ", "
                + to_string(this->nodes->size()));
        return;
    }

    double q_sum = 0;
    for (int i = 0; i < n_nodes; ++i) {
        const double scalar = this->get_scalar(i);
        if (this->nodes->get_marker(i) == TYPES.TETNODE)
            q_sum += scalar;
    }

    stringstream stream; stream << fixed << setprecision(3);
    stream << "Q / sum(" << this->scalar_label << ") = " << Q << " / " << q_sum << " = " << Q/q_sum;
    write_verbose_msg(stream.str());
}

// Calculate distance-dependent weights for a point with respect to the cell
template<int rank>
void TriangleInterpolator<rank>::get_weights(array<double,rank>& weights, const Point3 &point, const SimpleCell<rank>& scell) const {
    double w_sum = 0;
    for (int i = 0; i < rank; ++i) {
        double w = exp(decay_factor * point.distance(vertices[scell[i]]));
        weights[i] = w;
        w_sum += w;
    }

    // check for zero and nan
    require(w_sum > 0 && w_sum == w_sum, "Invalid interpolation weight: " + to_string(w_sum));
    w_sum = 1.0 / w_sum;
    for (int i = 0; i < rank; ++i)
        weights[i] *= w_sum;
}

// Interpolate both scalar and vector data inside or near the cell
template<int rank>
Solution TriangleInterpolator<rank>::interp_solution(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

    SimpleCell<rank> scell = cells[cell];

    // calculate weights or barycentric coordinates
    array<double,rank> weights;
    if (c >= 0) get_shape_functions(weights, Vec3(point), cell);
    else get_weights(weights, point, scell);

    // Interpolate electric field
    Vec3 vector_i(0.0);
    for (int i = 0; i < rank; ++i)
        vector_i += solutions[scell[i]].vector * weights[i];

    // Interpolate potential
    double scalar_i(0.0);
    for (int i = 0; i < rank; ++i)
        scalar_i += solutions[scell[i]].scalar * weights[i];

    return Solution(vector_i, scalar_i);
}

/* ==================================================================
 *  ====================== LinTriInterpolator ======================
 * ================================================================== */

void LinTriInterpolator::get_shape_functions(array<double,3>& sf, const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    const double v = tvec.dotProduct(pvec[face]);
    const double w = qvec.dotProduct(norms[face]);
    const double u = 1.0 - v - w;
    sf = {u, v, w};
}

SimpleCell<3> LinTriInterpolator::get_cell(const int tri) const {
    require(tri >= 0 && tri < faces->size(), "Invalid index: " + to_string(tri));
    return (*faces)[tri];
}

/* ==================================================================
 *  ====================== QuadTriInterpolator =====================
 * ================================================================== */

void QuadTriInterpolator::get_shape_functions(array<double,6>& sf, const Vec3& point, const int face) const {
    const Vec3 tvec = point - vert0[face];
    const Vec3 qvec = tvec.crossProduct(edge1[face]);
    const double v = tvec.dotProduct(pvec[face]);
    const double w = qvec.dotProduct(norms[face]);
    const double u = 1.0 - v - w;

    sf[0] = u * (2 * u - 1);
    sf[1] = v * (2 * v - 1);
    sf[2] = w * (2 * w - 1);
    sf[3] = 4 * u * v;
    sf[4] = 4 * v * w;
    sf[5] = 4 * w * u;

//    sf = {u, v, w, 0, 0, 0};
}

SimpleCell<6> QuadTriInterpolator::get_cell(const int tri) const {
    require(tri >= 0 && tri < faces->size(), "Invalid index: " + to_string(tri));
    if (mesh->quads.size() == 0)
        return QuadraticTri(0);

    const int n_quads_per_tri = 3;
    array<vector<unsigned>,3> edge_nodes;

    // locate hexahedral nodes that are located in the middle of edges
    for (int i = 0; i < n_quads_per_tri; ++i) {
        for (int quadnode : mesh->quads[n_quads_per_tri * tri + i])
            if (nodes->get_marker(quadnode) == TYPES.EDGECENTROID)
                edge_nodes[i].push_back(quadnode);
    }

    // find second order nodes
    const int n4 = common_entry(edge_nodes[0], edge_nodes[1]);
    const int n5 = common_entry(edge_nodes[1], edge_nodes[2]);
    const int n6 = common_entry(edge_nodes[2], edge_nodes[0]);

    return QuadraticTri((*faces)[tri], n4, n5, n6);
}

template class TriangleInterpolator<3> ;
template class TriangleInterpolator<6> ;
template class TetrahedronInterpolator<4> ;
template class TetrahedronInterpolator<10>;

} // namespace femocs
