/*
 * InterpolatorCells.cpp
 *
 *  Created on: 10.1.2018
 *      Author: veske
 */

#include "InterpolatorCells.h"
#include <float.h>

namespace femocs {

/* ==================================================================
 *  ====================== InterpolatorNodes =======================
 * ================================================================== */

InterpolatorNodes::InterpolatorNodes() :
        mesh(NULL), vec_label("vector"), norm_label("vector_norm"), scalar_label("scalar")
{
    reserve(0);
}

InterpolatorNodes::InterpolatorNodes(const string &vl, const string &nl, const string &sl) :
        mesh(NULL), vec_label(vl), norm_label(nl), scalar_label(sl)
{
    reserve(0);
}

void InterpolatorNodes::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));

    markers = vector<int>(N);
    solutions.clear();
    solutions.reserve(N);
}

void InterpolatorNodes::precompute() {
    const int n_nodes = mesh->nodes.size();
    reserve(n_nodes);

    for (int i = 0; i < n_nodes; ++i) {
        // Store node markers
        markers[i] = mesh->nodes.get_marker(i);
    }
}

void InterpolatorNodes::read(const string &file_name, const int flags) {
    string ftype = get_file_type(file_name);
    require(ftype == "restart", "Unimplemented file type: " + ftype);

    ifstream in(file_name);
    require(in.is_open(), "Can't open a file " + file_name);

    expect(flags, "No data will be read!");
    bool import_vec = flags & (1 << 2);
    bool import_norm = flags & (1 << 1);
    bool import_scalar = flags & (1 << 0);

    string label = "$";
    if (import_vec) label += vec_label;
    if (import_norm) label += norm_label;
    if (import_scalar) label += scalar_label;

    string str;
    while (in >> str) {
        if (str == label) {
            int n_nodes;
            in >> n_nodes >> GLOBALS.TIME >> GLOBALS.TIMESTEP;
            getline(in, str);

            reserve(n_nodes);

            Solution s(0);
            for (int ver = 0; ver < n_nodes; ++ver) {
                if (import_vec)
                    in.read(reinterpret_cast<char*>(&s.vector), sizeof(Vec3));
                if (import_norm)
                    in.read(reinterpret_cast<char*>(&s.norm), sizeof(double));
                if (import_scalar)
                    in.read(reinterpret_cast<char*>(&s.scalar), sizeof(double));
                append_solution(s);
            }
        }
    }

    in.close();
}

void InterpolatorNodes::write_xyz(ofstream& out) const {
    // write the beginning of xyz header
    FileWriter::write_xyz(out);

    // write the header for Ovito
    out << "properties=id:I:1:pos:R:3:marker:I:1:" <<
            "force:R:3:" << norm_label << ":R:1:" << scalar_label << ":R:1" << endl;

    // write data
    const int n_nodes = size();
    for (int i = 0; i < n_nodes; ++i)
        out << i << " " << get_vertex(i) << " " << markers[i] << " " << solutions[i] << endl;
}

void InterpolatorNodes::write_bin(ofstream &out) const {
    for (Solution const &s : solutions)
            out.write ((char*)&s.scalar, sizeof (double));
}

void InterpolatorNodes::write_vtk_points_and_cells(ofstream& out) const {
    const int n_nodes = size();
    const int n_cells = size();
    const int dim = 1;

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << get_vertex(i) << "\n";

    // Output # vertices and vertex indices
    out << "CELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << i << "\n";

    // Output cell types
    out << "CELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << TYPES.VTK.VERTEX << "\n";
}

void InterpolatorNodes::write_vtk_point_data(ofstream& out) const {
    const int n_nodes = size();
    out << "POINT_DATA " << n_nodes << "\n";

    // write node IDs
    out << "SCALARS node-ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << i << "\n";

    // write node markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << markers[i] << "\n";

    // write vector norm data
    out << "SCALARS " + norm_label + " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << solutions[i].norm << "\n";

    // write scalar data
    out << "SCALARS " + scalar_label + " double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i)
        out << solutions[i].scalar << "\n";

    // write vector data
    out << "VECTORS " << vec_label << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << solutions[i].vector << "\n";
}

void InterpolatorNodes::print_statistics() const {
    if (!MODES.VERBOSE) return;

    const int n_atoms = size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        if ( get_vector_norm(i) == 0) continue;
        double s = get_scalar(i);
        Vec3 v = get_vector(i);

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

double InterpolatorNodes::max_norm() const {
    double max_norm = -DBL_MAX;
    for (Solution s : solutions)
        max_norm = max(max_norm, s.norm);
    return max_norm;
}

/* ==================================================================
 *  ====================== InterpolatorCells =======================
 * ================================================================== */

template<int dim>
void InterpolatorCells<dim>::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));

    centroids.clear();
    centroids.reserve(N);
    markers = vector<int>(N);
    neighbours = vector<vector<int>>(N);
}

template<int dim>
void InterpolatorCells<dim>::write_vtk_points_and_cells(ofstream& out) const {
    const int n_nodes = nodes->size();
    const int celltype = get_cell_type();
    const int n_cells = size();

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << nodes->get_vertex(i) << "\n";

    // Output # vertices and vertex indices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << get_cell(i) << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";
}

template<int dim>
void InterpolatorCells<dim>::write_vtk_point_data(ofstream &out) const {
    nodes->write_vtk_point_data(out);
}

template<int dim>
void InterpolatorCells<dim>::write_vtk_cell_data(ofstream& out) const {
    const int n_cells = size();
    out << "CELL_DATA " << n_cells << "\n";

    // write cell IDs
    out << "SCALARS cell-ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_cells; ++i)
        out << i << "\n";

    // write cell markers
    out << "SCALARS cell-marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_cells; ++i)
        out << markers[i] << "\n";
}

template<int dim>
int InterpolatorCells<dim>::locate_cell(const Point3 &point, const int cell_guess) const {
    const int n_cells = neighbours.size();
    require(cell_guess < n_cells, "Index out of bounds: " + d2s(cell_guess));

    if (cell_guess >= 0) {
        // === Check the guessed cell
        if (point_in_cell(point, cell_guess))
            return cell_guess;

        // === Check if point is surrounded by one of the neighbouring cells
        for (int cell : neighbours[cell_guess]) {
            if (point_in_cell(point, cell))
                return cell;
        }
    }

    // === In case of no success, loop through all the cells
    double min_distance2 = 1e100;
    int min_index = 0;

    for (int cell = 0; cell < n_cells; ++cell) {
        // If correct cell is found, we're done
        if (markers[cell] == 0 && point_in_cell(point, cell))
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

template<int dim>
Solution InterpolatorCells<dim>::interp_solution(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < size(), "Index out of bounds: " + d2s(cell));

    // calculate interpolation weights (shape function of dummy nodal distance based weights)
    array<double,dim> weights = shape_functions(Vec3(point), cell);

    SimpleCell<dim> scell = get_cell(cell);

    // Interpolate vector data
    Vec3 vector_i(0.0);
    for (int i = 0; i < dim; ++i)
        vector_i += nodes->get_vector(scell[i]) * weights[i];

    // Interpolate scalar data
    double scalar_i(0.0);
    for (int i = 0; i < dim; ++i)
        scalar_i += nodes->get_scalar(scell[i]) * weights[i];

    return Solution(vector_i, scalar_i);
}

template<int dim>
Solution InterpolatorCells<dim>::interp_solution_v2(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < size(), "Index out of bounds: " + d2s(cell));

    // calculate shape functions and their gradients
    array<double,dim> sf = shape_functions(point, cell);
    array<Vec3,dim> sfg = shape_fun_grads(point, cell);

    // using them as weights, interpolate scalar data and its gradient
    SimpleCell<dim> scell = get_cell(cell);
    Vec3 vector_i(0.0);
    double scalar_i(0.0);

    for (int i = 0; i < dim; ++i) {
        double scalar = nodes->get_scalar(scell[i]);
        vector_i -= sfg[i] * scalar;
        scalar_i += sf[i] * scalar;
    }

    return Solution(vector_i, scalar_i);
}

template<int dim>
Vec3 InterpolatorCells<dim>::interp_gradient(const Point3 &point, const int cell) const {
    require(cell >= 0 && cell < size(), "Index out of bounds: " + d2s(cell));

    // calculate shape function gradients
    array<Vec3, dim> sfg = shape_fun_grads(point, cell);

    // using them as weights, interpolate the gradient of scalar data
    SimpleCell<dim> scell = get_cell(cell);
    Vec3 vector_i(0.0);

    for (int i = 0; i < dim; ++i)
        vector_i -= sfg[i] * nodes->get_scalar(scell[i]);

    return vector_i;
}

template<int dim>
Vec3 InterpolatorCells<dim>::interp_gradient(const int cell, const int node) const {
    require(cell >= 0 && cell < size(), "Index out of bounds: " + d2s(cell));

    // calculate shape function gradients
    array<Vec3,dim> sfg = shape_fun_grads(cell, node);

    // using them as weights, interpolate the gradient of scalar data
    SimpleCell<dim> scell = get_cell(cell);
    Vec3 vector_i(0.0);

    for (int i = 0; i < dim; ++i)
        vector_i -= sfg[i] * nodes->get_scalar(scell[i]);

    return vector_i;
}

template<int dim>
Solution InterpolatorCells<dim>::locate_interpolate(const Point3 &point, int& cell) const {
    cell = locate_cell(point, abs(cell));
    return interp_solution(point, cell);
}

template<int dim>
Solution InterpolatorCells<dim>::locate_interpolate_v2(const Point3 &point, int& cell) const {
    cell = locate_cell(point, abs(cell));
    return interp_solution_v2(point, cell);
}

template<int dim>
int InterpolatorCells<dim>::common_entry(vector<unsigned>& vec1, vector<unsigned>& vec2) const {
    for (unsigned i : vec1)
        for (unsigned j : vec2)
            if (i == j) return i;
    return -1;
}

template<int dim>
double InterpolatorCells<dim>::determinant(const Vec3 &v1, const Vec3 &v2) const {
    return v1.x * (v2.y - v2.z) - v1.y * (v2.x - v2.z) + v1.z * (v2.x - v2.y);
}

template<int dim>
double InterpolatorCells<dim>::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const {
    return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z)
            + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

template<int dim>
double InterpolatorCells<dim>::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) const {
    const double det1 = determinant(v2, v3, v4);
    const double det2 = determinant(v1, v3, v4);
    const double det3 = determinant(v1, v2, v4);
    const double det4 = determinant(v1, v2, v3);

    return det4 - det3 + det2 - det1;
}

template<int dim>
double InterpolatorCells<dim>::determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) const {
    double det1 = determinant(Vec3(v1.y,v1.z,v1.w), Vec3(v2.y,v2.z,v2.w), Vec3(v3.y,v3.z,v3.w));
    double det2 = determinant(Vec3(v1.x,v1.z,v1.w), Vec3(v2.x,v2.z,v2.w), Vec3(v3.x,v3.z,v3.w));
    double det3 = determinant(Vec3(v1.x,v1.y,v1.w), Vec3(v2.x,v2.y,v2.w), Vec3(v3.x,v3.y,v3.w));
    double det4 = determinant(Vec3(v1.x,v1.y,v1.z), Vec3(v2.x,v2.y,v2.z), Vec3(v3.x,v3.y,v3.z));

    return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
}

/* ==================================================================
 *  ====================== LinearTetrahedra ========================
 * ================================================================== */

LinearTetrahedra::LinearTetrahedra() :
        InterpolatorCells<4>(), tets(NULL) {}

LinearTetrahedra::LinearTetrahedra(const InterpolatorNodes* n) :
        InterpolatorCells<4>(n), tets(NULL) {}

void LinearTetrahedra::reserve(const int N) {
    InterpolatorCells<4>::reserve(N);

    det0.clear(); det0.reserve(N);
    det1.clear(); det1.reserve(N);
    det2.clear(); det2.reserve(N);
    det3.clear(); det3.reserve(N);
    det4.clear(); det4.reserve(N);
    tet_not_valid.clear(); tet_not_valid.reserve(N);
}

void LinearTetrahedra::precompute() {
    require(mesh && tets, "NULL pointers can't be used!");
    const int n_elems = tets->size();
    const int n_nodes = mesh->nodes.size();

    expect(n_nodes > 0 && n_elems > 0, "Interpolator expects non-empty mesh!");
    double d0, d1, d2, d3, d4;

    reserve(n_elems);

    // Store the constant for smoothing
    decay_factor = -1.0 / tets->stat.edgemax;

    // Calculate which tetrahedra are connected to the node
    vector<vector<int>> node2tets(n_nodes);
    for (int tet = 0; tet < n_elems; ++tet) {
        for (int node : (*tets)[tet])
            node2tets[node].push_back(tet);
    }

    // loop through all the tetrahedra
    for (int tet = 0; tet < n_elems; ++tet) {
        SimpleElement selem = (*tets)[tet];

        // store nearest neighbours of tetrahedron
        vector<int> nnbors = tets->get_neighbours(tet);
        for (int nbor : nnbors)
            if (nbor >= 0)
                neighbours[tet].push_back(nbor);

        // store next nearest neighbours of tetrahedron
        for (int node : selem)
            for (int nbor_tet : node2tets[node]) {
                if (nbor_tet != tet && nbor_tet != nnbors[0] && nbor_tet != nnbors[1]
                         && nbor_tet != nnbors[2] && nbor_tet != nnbors[3])
                    neighbours[tet].push_back(nbor_tet);
            }

        // Calculate centroids of tetrahedra
        centroids.push_back(tets->get_centroid(tet));

        /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
         * barycentric coordinate of tetrahedra using the relations below */

        Vec3 v1 = mesh->nodes.get_vec(selem[0]);
        Vec3 v2 = mesh->nodes.get_vec(selem[1]);
        Vec3 v3 = mesh->nodes.get_vec(selem[2]);
        Vec3 v4 = mesh->nodes.get_vec(selem[3]);

        /* =====================================================================================
         * det0 = |x1 y1 z1 1|
                  |x2 y2 z2 1|
                  |x3 y3 z3 1|
                  |x4 y4 z4 1|  */
        d0 = determinant(v1, v2, v3, v4);
        tet_not_valid.push_back(fabs(d0) < zero);
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

bool LinearTetrahedra::point_in_cell(const Vec3& point, const int i) const {
    require(i >= 0 && i < (int)det0.size(), "Index out of bounds: " + d2s(i));

    // Ignore co-planar tetrahedra
    // no need to check because Tetgen guarantees non-co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    const Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < -zero) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < -zero) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < -zero) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < -zero) return false;

    // All bcc-s are >= 0, so point is inside the tetrahedron
    return true;
}

array<double,4> LinearTetrahedra::shape_functions(const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < (int)det0.size(), "Index out of bounds: " + d2s(tet));

    const Vec4 pt(point, 1);
    return {
        zero + det0[tet] * pt.dotProduct(det1[tet]),
        zero + det0[tet] * pt.dotProduct(det2[tet]),
        zero + det0[tet] * pt.dotProduct(det3[tet]),
        zero + det0[tet] * pt.dotProduct(det4[tet])
    };
}

array<Vec3,4> LinearTetrahedra::shape_fun_grads(const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < size(), "Index out of bounds: " + d2s(tet));

    // The gradient of shape function inside linear tetrahedron does not depend
    // on the point of interest, therefore its value could be precalculated and
    // used upon request. However, for the consistency with other cells it's not
    // done now. Maybe in the future, if someone wants to use this function very heavily.

    // calculate 4x4 Jacobian
    // the upmost row consists of 1-s and is not used, therefore no need to store it
    SimpleCell<4> cell = get_cell(tet);
    array<Vec3, 4> J = {
        mesh->nodes[cell[0]],
        mesh->nodes[cell[1]],
        mesh->nodes[cell[2]],
        mesh->nodes[cell[3]]
    };

    // calculate determinant (and its inverse) of Jacobian
    double Jdet = -1.0 * determinant(J[0], J[1], J[2], J[3]);
    require(fabs(Jdet) > 1e-15, "Singular Jacobian can't be handled!");
    Jdet = 1.0 / Jdet;

    // calculate inverse of Jacobian
    // The first column is omitted, as it's not used
    Vec3 J01 = J[0] - J[1];
    Vec3 J02 = J[0] - J[2];
    Vec3 J03 = J[0] - J[3];
    Vec3 J12 = J[1] - J[2];
    Vec3 J13 = J[1] - J[3];
    Vec3 J23 = J[2] - J[3];

    array<Vec3, 4> Jinv =
    { Vec3(
          J[2].y*J13.z - J[3].y*J12.z - J[1].y*J23.z,
         -J[2].x*J13.z + J[3].x*J12.z + J[1].x*J23.z,
          J[2].x*J13.y - J[3].x*J12.y - J[1].x*J23.y
    ), Vec3(
         -J[2].y*J03.z + J[3].y*J02.z + J[0].y*J23.z,
          J[2].x*J03.z - J[3].x*J02.z - J[0].x*J23.z,
         -J[2].x*J03.y + J[3].x*J02.y + J[0].x*J23.y
    ), Vec3(
          J[1].y*J03.z - J[3].y*J01.z - J[0].y*J13.z,
         -J[1].x*J03.z + J[3].x*J01.z + J[0].x*J13.z,
          J[1].x*J03.y - J[3].x*J01.y - J[0].x*J13.y
    ), Vec3(
         -J[1].y*J02.z + J[2].y*J01.z + J[0].y*J12.z,
          J[1].x*J02.z - J[2].x*J01.z - J[0].x*J12.z,
         -J[1].x*J02.y + J[2].x*J01.y + J[0].x*J12.y
    )};
    for (int i = 0; i < 4; ++i)
        Jinv[i] *= Jdet;

    // calculate gradient of shape functions in xyz-space
    return {Jinv[0], Jinv[1], Jinv[2], Jinv[3]};
}

void LinearTetrahedra::narrow_search_to(const int region) {
    const int n_cells = size();
    require(mesh->tets.size() == n_cells,
            "Mismatch between tetrahedral mesh and interpolator sizes (" + d2s(mesh->tets.size()) + " vs " + d2s(n_cells)
            + ")\nindicates, that LinearTetrahedra are not properly pre-computed!");

    const int surf_end = mesh->nodes.indxs.surf_end;

    if (region == TYPES.VACUUM) {
        for (int i = 0; i < n_cells; ++i)
            markers[i] = tets->get_marker(i) != TYPES.VACUUM;
    } else if (region == TYPES.BULK) {
        for (int i = 0; i < n_cells; ++i)
            markers[i] = tets->get_marker(i) == TYPES.VACUUM;
    } else if (region == TYPES.SURFACE) {
        for (int i = 0; i < n_cells; ++i) {
            markers[i] = get_cell(i) > surf_end;
        }
    } else if (region == TYPES.NONE) {
        markers = vector<int>(n_cells);
    } else
        require(false, "Unimplemented region: " + d2s(region));
}

/* ==================================================================
 *  ===================== QuadraticTetrahedra ======================
 * ================================================================== */

QuadraticTetrahedra::QuadraticTetrahedra() :
        InterpolatorCells<10>(), tets(NULL), lintet(NULL) {}

QuadraticTetrahedra::QuadraticTetrahedra(const InterpolatorNodes* n, const LinearTetrahedra* l) :
    InterpolatorCells<10>(n), tets(NULL), lintet(l) {}

void QuadraticTetrahedra::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));
    markers = vector<int>(N);
    cells.clear();
    cells.reserve(N);
}

void QuadraticTetrahedra::precompute() {
    require(mesh && tets && lintet, "NULL pointers can't be used!");
    require(tets->size() == lintet->size(),
            "Mismatch between tetrahedral mesh and interpolator sizes (" + d2s(tets->size()) + " vs " + d2s(lintet->size())
            + ")\nindicates, that LinearTetrahedra are not properly pre-computed!");

    const int n_elems = tets->size();
    reserve(n_elems);

    // Store the constant for smoothing
    this->decay_factor = -1.0 / tets->stat.edgemax;

    // Loop through all the tetrahedra
    for (int i = 0; i < n_elems; ++i) {
        // make the markers to correspond to tetrahedra
        markers[i] = tets->get_marker(i);

        // Calculate and store 10-noded tetrahedra
        cells.push_back(calc_cell(i));
    }
}

bool QuadraticTetrahedra::point_in_cell(const Vec3& point, const int cell) const {
    return lintet->point_in_cell(point, cell);
}

int QuadraticTetrahedra::locate_cell(const Point3 &point, const int cell_guess) const {
    return lintet->locate_cell(point, cell_guess);
}

array<double,10> QuadraticTetrahedra::shape_functions(const Vec3& point, const int tet) const {
    array<double,4> bcc = lintet->shape_functions(point, tet);

    const double b1 = bcc[0];
    const double b2 = bcc[1];
    const double b3 = bcc[2];
    const double b4 = bcc[3];

    return {
        b1 * (2 * b1 - 1),
        b2 * (2 * b2 - 1),
        b3 * (2 * b3 - 1),
        b4 * (2 * b4 - 1),
        4 * b1 * b2,
        4 * b2 * b3,
        4 * b3 * b1,
        4 * b1 * b4,
        4 * b2 * b4,
        4 * b3 * b4
    };
//    return {b1, b2, b3, b4, 0, 0, 0, 0, 0, 0};
}

/*
 * Calculate gradient of shape function for 10-noded tetrahedra.
 * The theory can be found from lecture notes at
 * https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
 *
 * The summary:
 * a) As the shape functions of tetrahedra are defined by 4 barycentric coordinates (bcc),
 * the coordinate matrix must be appended with additional row of 1-s:
 *
 *    double xyz[4][10] = {{1 x0 y0 z0}, {1 x1 y1 z1}, ... {1 x9 y9 z9}},
 *
 * b) Calculate gradient of shape functions in bcc-space:
 *    double dN[4][10] = {
 *        {4*bcc[0]-1, 0, 0, 0,  4*bcc[1],        0, 4*bcc[2], 4*bcc[3],        0,        0},
 *        {0, 4*bcc[1]-1, 0, 0,  4*bcc[0], 4*bcc[2],        0,        0, 4*bcc[3],        0},
 *        {0, 0, 4*bcc[2]-1, 0,         0, 4*bcc[1], 4*bcc[0],        0,        0, 4*bcc[3]},
 *        {0, 0, 0, 4*bcc[3]-1,         0,        0,        0, 4*bcc[0], 4*bcc[1], 4*bcc[2]}
 *    };
 *
 * c) Calculate 4x4 Jacobian, upmost row of J will consist of 1-s
 *    double J[4][4] = {0};
 *    for (int k = 0; k < 10; ++k)
 *      for (int i = 0; i < 4; ++i)
 *        for (int j = 0; j < 4; ++j)
 *          J[i][j] += dN[i][k] * xyz[j][k];
 *
 * d) Calculate determinant and inverse of Jacobian in a standard way;
 *    omit leftmost column of J inverse, as it is not needed; obtain array<Vec3,4> Jinv.
 *
 * e) Calculate gradient of shape functions in xyz-space
 *    array<Vec3, 10> sfg;
 *    for (int k = 0; k < 10; ++k)
 *      for (int j = 0; j < 4; ++j)
 *        sfg[k] += Jinv[j] * dN[j][k];
 *
 * As dN contains many zeros and the upmost row of J and leftmost column of Jinv is not used,
 * the above calculations can be optimized to reduce nr of needed flops.
 * For instance there's no need to calculate dN explicitly.
 */
array<Vec3,10> QuadraticTetrahedra::shape_fun_grads(const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < (int)cells.size(), "Index out of bounds: " + d2s(tet));

    array<double,4> bcc = lintet->shape_functions(point, tet);

    // Store tetrahedron nodes for more convenient access
    SimpleCell<10> cell = get_cell(tet);
    array<Vec3, 10> node;
    for (int i = 0; i < 10; ++i)
        node[i] = mesh->nodes[cell[i]];

    // calculate 4x4 Jacobian
    // the upmost row consists of 1-s and is not used, therefore no need to store it
    array<Vec3, 4> J = {
        node[0]*(bcc[0]-0.25) + node[4]*bcc[1] + node[6]*bcc[2] + node[7]*bcc[3],
        node[4]*bcc[0] + node[1]*(bcc[1]-0.25) + node[5]*bcc[2] + node[8]*bcc[3],
        node[6]*bcc[0] + node[5]*bcc[1] + node[2]*(bcc[2]-0.25) + node[9]*bcc[3],
        node[7]*bcc[0] + node[8]*bcc[1] + node[9]*bcc[2] + node[3]*(bcc[3]-0.25)
    };
    for (int i = 0; i < 4; ++i)
        J[i] *= 4.0;

    // calculate determinant (and its inverse) of Jacobian
    double Jdet = -1.0 * determinant(J[0], J[1], J[2], J[3]);
    require(fabs(Jdet) > 1e-15, "Singular Jacobian can't be handled!");
    Jdet = 1.0 / Jdet;

    // calculate inverse of Jacobian
    // The first column is omitted, as it's not used
    Vec3 J01 = J[0] - J[1];
    Vec3 J02 = J[0] - J[2];
    Vec3 J03 = J[0] - J[3];
    Vec3 J12 = J[1] - J[2];
    Vec3 J13 = J[1] - J[3];
    Vec3 J23 = J[2] - J[3];

    array<Vec3, 4> Jinv =
    { Vec3(
          J[2].y*J13.z - J[3].y*J12.z - J[1].y*J23.z,
         -J[2].x*J13.z + J[3].x*J12.z + J[1].x*J23.z,
          J[2].x*J13.y - J[3].x*J12.y - J[1].x*J23.y
    ), Vec3(
         -J[2].y*J03.z + J[3].y*J02.z + J[0].y*J23.z,
          J[2].x*J03.z - J[3].x*J02.z - J[0].x*J23.z,
         -J[2].x*J03.y + J[3].x*J02.y + J[0].x*J23.y
    ), Vec3(
          J[1].y*J03.z - J[3].y*J01.z - J[0].y*J13.z,
         -J[1].x*J03.z + J[3].x*J01.z + J[0].x*J13.z,
          J[1].x*J03.y - J[3].x*J01.y - J[0].x*J13.y
    ), Vec3(
         -J[1].y*J02.z + J[2].y*J01.z + J[0].y*J12.z,
          J[1].x*J02.z - J[2].x*J01.z - J[0].x*J12.z,
         -J[1].x*J02.y + J[2].x*J01.y + J[0].x*J12.y
    )};
    for (int i = 0; i < 4; ++i)
        Jinv[i] *= Jdet;

    // calculate gradient of shape functions in xyz-space
    return {
        Jinv[0] * (4*bcc[0]-1),
        Jinv[1] * (4*bcc[1]-1),
        Jinv[2] * (4*bcc[2]-1),
        Jinv[3] * (4*bcc[3]-1),
        (Jinv[0]*bcc[1] + Jinv[1]*bcc[0]) * 4.0,
        (Jinv[1]*bcc[2] + Jinv[2]*bcc[1]) * 4.0,
        (Jinv[0]*bcc[2] + Jinv[2]*bcc[0]) * 4.0,
        (Jinv[0]*bcc[3] + Jinv[3]*bcc[0]) * 4.0,
        (Jinv[1]*bcc[3] + Jinv[3]*bcc[1]) * 4.0,
        (Jinv[2]*bcc[3] + Jinv[3]*bcc[2]) * 4.0
    };
}

array<Vec3,10> QuadraticTetrahedra::shape_fun_grads_slow(const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < (int)cells.size(), "Index out of bounds: " + d2s(tet));

    array<double,4> bcc = lintet->shape_functions(point, tet);

    cout << fixed << setprecision(3);

    cout << "bcc = ";
    for (int i = 0; i < 4; ++i)
        cout << bcc[i] << " ";
    cout << endl;

    // Store tetrahedron nodes for more convenient access
    SimpleCell<10> cell = cells[tet];

    array<Vec4,10> xyz;
    for (int i = 0; i < 10; ++i) {
        Vec3 node = mesh->nodes[cell[i]];
        xyz[i][0] = 1;
        for (int j = 0; j < 3; ++j)
            xyz[i][j+1] = node[j];
    }

    cout << "xyz = \n";
    for (int i = 0; i < 10; ++i)
        cout << xyz[i] << endl;

    // Calculate gradient of shape functions in bcc-space:
    array<Vec4,10> dN = {
            Vec4(4*bcc[0]-1, 0, 0, 0),
            Vec4(0, 4*bcc[1]-1, 0, 0),
            Vec4(0, 0, 4*bcc[2]-1, 0),
            Vec4(0, 0, 0, 4*bcc[3]-1),
            Vec4(4*bcc[1], 4*bcc[0], 0, 0),
            Vec4(0, 4*bcc[2], 4*bcc[1], 0),
            Vec4(4*bcc[2], 0, 4*bcc[0], 0),
            Vec4(4*bcc[3], 0, 0, 4*bcc[0]),
            Vec4(0, 4*bcc[3], 0, 4*bcc[1]),
            Vec4(0, 0, 4*bcc[3],4*bcc[2])
    };

    cout << "dN = \n";
    for (int i = 0; i < 10; ++i)
        cout << dN[i] << endl;

    // calculate 4x4 Jacobian
    array<Vec4,4> J;
    for (int k = 0; k < 10; ++k)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                J[j][i] += dN[k][i] * xyz[k][j];

    cout << "J = {\n";
    for (int i = 0; i < 4; ++i) {
        cout << "{";
        for (int j = 0; j < 3; ++j)
            cout << J[i][j] << ", ";
        if (i != 3)
            cout << J[i][3] << "}," << endl;
        else
            cout << J[i][3] << "}\n}" << endl;
    }

    // calculate determinant (and its inverse) of Jacobian
    double Jdet = determinant(J[0], J[1], J[2], J[3]);

    cout << "Jdet = " << Jdet << endl;

    require(fabs(Jdet) > 1e-15, "Singular Jacobian can't be handled!");
    Jdet = 1.0 / Jdet;

    array<Vec4, 4> Jinv;
    Jinv[0][0] = J[1][1]  * J[2][2] * J[3][3] -
             J[1][1]  * J[2][3] * J[3][2] -
             J[2][1]  * J[1][2]  * J[3][3] +
             J[2][1]  * J[1][3]  * J[3][2] +
             J[3][1] * J[1][2]  * J[2][3] -
             J[3][1] * J[1][3]  * J[2][2];

    Jinv[1][0] = -J[1][0]  * J[2][2] * J[3][3] +
              J[1][0]  * J[2][3] * J[3][2] +
              J[2][0]  * J[1][2]  * J[3][3] -
              J[2][0]  * J[1][3]  * J[3][2] -
              J[3][0] * J[1][2]  * J[2][3] +
              J[3][0] * J[1][3]  * J[2][2];

    Jinv[2][0] = J[1][0]  * J[2][1] * J[3][3] -
             J[1][0]  * J[2][3] * J[3][1] -
             J[2][0]  * J[1][1] * J[3][3] +
             J[2][0]  * J[1][3] * J[3][1] +
             J[3][0] * J[1][1] * J[2][3] -
             J[3][0] * J[1][3] * J[2][1];

    Jinv[3][0] = -J[1][0]  * J[2][1] * J[3][2] +
               J[1][0]  * J[2][2] * J[3][1] +
               J[2][0]  * J[1][1] * J[3][2] -
               J[2][0]  * J[1][2] * J[3][1] -
               J[3][0] * J[1][1] * J[2][2] +
               J[3][0] * J[1][2] * J[2][1];

    Jinv[0][1] = -J[0][1]  * J[2][2] * J[3][3] +
              J[0][1]  * J[2][3] * J[3][2] +
              J[2][1]  * J[0][2] * J[3][3] -
              J[2][1]  * J[0][3] * J[3][2] -
              J[3][1] * J[0][2] * J[2][3] +
              J[3][1] * J[0][3] * J[2][2];

    Jinv[1][1] = J[0][0]  * J[2][2] * J[3][3] -
             J[0][0]  * J[2][3] * J[3][2] -
             J[2][0]  * J[0][2] * J[3][3] +
             J[2][0]  * J[0][3] * J[3][2] +
             J[3][0] * J[0][2] * J[2][3] -
             J[3][0] * J[0][3] * J[2][2];

    Jinv[2][1] = -J[0][0]  * J[2][1] * J[3][3] +
              J[0][0]  * J[2][3] * J[3][1] +
              J[2][0]  * J[0][1] * J[3][3] -
              J[2][0]  * J[0][3] * J[3][1] -
              J[3][0] * J[0][1] * J[2][3] +
              J[3][0] * J[0][3] * J[2][1];

    Jinv[3][1] = J[0][0]  * J[2][1] * J[3][2] -
              J[0][0]  * J[2][2] * J[3][1] -
              J[2][0]  * J[0][1] * J[3][2] +
              J[2][0]  * J[0][2] * J[3][1] +
              J[3][0] * J[0][1] * J[2][2] -
              J[3][0] * J[0][2] * J[2][1];

    Jinv[0][2] = J[0][1]  * J[1][2] * J[3][3] -
             J[0][1]  * J[1][3] * J[3][2] -
             J[1][1]  * J[0][2] * J[3][3] +
             J[1][1]  * J[0][3] * J[3][2] +
             J[3][1] * J[0][2] * J[1][3] -
             J[3][1] * J[0][3] * J[1][2];

    Jinv[1][2] = -J[0][0]  * J[1][2] * J[3][3] +
              J[0][0]  * J[1][3] * J[3][2] +
              J[1][0]  * J[0][2] * J[3][3] -
              J[1][0]  * J[0][3] * J[3][2] -
              J[3][0] * J[0][2] * J[1][3] +
              J[3][0] * J[0][3] * J[1][2];

    Jinv[2][2] = J[0][0]  * J[1][1] * J[3][3] -
              J[0][0]  * J[1][3] * J[3][1] -
              J[1][0]  * J[0][1] * J[3][3] +
              J[1][0]  * J[0][3] * J[3][1] +
              J[3][0] * J[0][1] * J[1][3] -
              J[3][0] * J[0][3] * J[1][1];

    Jinv[3][2] = -J[0][0]  * J[1][1] * J[3][2] +
               J[0][0]  * J[1][2] * J[3][1] +
               J[1][0]  * J[0][1] * J[3][2] -
               J[1][0]  * J[0][2] * J[3][1] -
               J[3][0] * J[0][1] * J[1][2] +
               J[3][0] * J[0][2] * J[1][1];

    Jinv[0][3] = -J[0][1] * J[1][2] * J[2][3] +
              J[0][1] * J[1][3] * J[2][2] +
              J[1][1] * J[0][2] * J[2][3] -
              J[1][1] * J[0][3] * J[2][2] -
              J[2][1] * J[0][2] * J[1][3] +
              J[2][1] * J[0][3] * J[1][2];

    Jinv[1][3] = J[0][0] * J[1][2] * J[2][3] -
             J[0][0] * J[1][3] * J[2][2] -
             J[1][0] * J[0][2] * J[2][3] +
             J[1][0] * J[0][3] * J[2][2] +
             J[2][0] * J[0][2] * J[1][3] -
             J[2][0] * J[0][3] * J[1][2];

    Jinv[2][3] = -J[0][0] * J[1][1] * J[2][3] +
               J[0][0] * J[1][3] * J[2][1] +
               J[1][0] * J[0][1] * J[2][3] -
               J[1][0] * J[0][3] * J[2][1] -
               J[2][0] * J[0][1] * J[1][3] +
               J[2][0] * J[0][3] * J[1][1];

    Jinv[3][3] = J[0][0] * J[1][1] * J[2][2] -
              J[0][0] * J[1][2] * J[2][1] -
              J[1][0] * J[0][1] * J[2][2] +
              J[1][0] * J[0][2] * J[2][1] +
              J[2][0] * J[0][1] * J[1][2] -
              J[2][0] * J[0][2] * J[1][1];
    for (int i = 0; i < 4; ++i)
        Jinv[i] *= Jdet;

    cout << "Jinv = \n";
    for (int i = 0; i < 4; ++i)
        cout << Jinv[i] << endl;

    array<Vec3,10> sfg;

    // calculate gradient of shape functions in xyz-space
    for (int k = 0; k < 10; ++k) {
        sfg[k] = Vec3(0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j)
                sfg[k][i] += Jinv[j][i+1] * dN[k][j];
    }

    return sfg;
}

void QuadraticTetrahedra::test_shape_funs() {
    array<Vec3, 10> sfg1, sfg2;

    cout << fixed << setprecision(3);

    int tet = 0;
    QuadraticTet cell = cells[tet];
    for (int i = 0; i < 4; ++i) {
        cout << "\ni = " << i << endl;
        Point3 point = mesh->nodes[cell[i]];
        sfg1 = shape_fun_grads(point, tet);
        sfg2 = shape_fun_grads_slow(point, tet);

        cout << "delta sfg = " << endl;
        for (int i = 0; i < 10; ++i)
            cout << sfg1[i] - sfg2[i] << endl;
    }
}

SimpleCell<10> QuadraticTetrahedra::calc_cell(const int tet) const {
    if (mesh->hexs.size() <= tet)
        return QuadraticTet(0);

    const int n_hexs_per_tet = 4;
    array<vector<unsigned>,4> edge_nodes;

    // locate hexahedral nodes that are located in the middle of edges
    for (int i = 0; i < n_hexs_per_tet; ++i) {
        for (int hexnode : mesh->hexs[n_hexs_per_tet * tet + i])
            if (mesh->nodes.get_marker(hexnode) == TYPES.EDGECENTROID)
                edge_nodes[i].push_back(hexnode);
    }

    // find second order nodes
    const int n5 = common_entry(edge_nodes[0], edge_nodes[1]);
    const int n6 = common_entry(edge_nodes[1], edge_nodes[2]);
    const int n7 = common_entry(edge_nodes[2], edge_nodes[0]);
    const int n8 = common_entry(edge_nodes[0], edge_nodes[3]);
    const int n9 = common_entry(edge_nodes[1], edge_nodes[3]);
    const int n10= common_entry(edge_nodes[2], edge_nodes[3]);

    return QuadraticTet((*tets)[tet], n5, n6, n7, n8, n9, n10);
}

/* ==================================================================
 *  ======================= LinearHexahedra ========================
 * ================================================================== */

LinearHexahedra::LinearHexahedra() :
        InterpolatorCells<8>(), hexs(NULL), lintet(NULL) {}

LinearHexahedra::LinearHexahedra(const InterpolatorNodes* n, const LinearTetrahedra* l) :
        InterpolatorCells<8>(n), hexs(NULL), lintet(l) {}

void LinearHexahedra::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));

    markers = vector<int>(N);

    f0s.clear(); f0s.reserve(N);
    f1s.clear(); f1s.reserve(N);
    f2s.clear(); f2s.reserve(N);
    f3s.clear(); f3s.reserve(N);
    f4s.clear(); f4s.reserve(N);
    f5s.clear(); f5s.reserve(N);
    f6s.clear(); f6s.reserve(N);
    f7s.clear(); f7s.reserve(N);
}

void LinearHexahedra::precompute() {
    require(mesh && hexs && lintet, "NULL pointers can't be used!");
    require(mesh->tets.size() == lintet->size(),
            "Mismatch between tetrahedral mesh and interpolator sizes (" + d2s(mesh->tets.size()) + " vs " + d2s(lintet->size())
            + ")\nindicates, that LinearTetrahedra are not properly pre-computed!");

    const int n_elems = hexs->size();
    const int n_nodes = mesh->nodes.size();
    expect(n_elems > 0 && n_nodes > 0, "Interpolator expects non-empty mesh!");
    reserve(n_elems);

    decay_factor = -1.0 / mesh->tets.stat.edgemax;

    // Loop through all the hexahedra
    for (int hex = 0; hex < n_elems; ++hex) {
        // make the markers to correspond to hexahedra
        markers[hex] = hexs->get_marker(hex);

        // pre-calculate data to make iterpolation faster
        SimpleHex shex = (*hexs)[hex];
        const Vec3 x1 = mesh->nodes.get_vec(shex[0]);
        const Vec3 x2 = mesh->nodes.get_vec(shex[1]);
        const Vec3 x3 = mesh->nodes.get_vec(shex[2]);
        const Vec3 x4 = mesh->nodes.get_vec(shex[3]);
        const Vec3 x5 = mesh->nodes.get_vec(shex[4]);
        const Vec3 x6 = mesh->nodes.get_vec(shex[5]);
        const Vec3 x7 = mesh->nodes.get_vec(shex[6]);
        const Vec3 x8 = mesh->nodes.get_vec(shex[7]);

        f0s.push_back( Vec3((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8) / 8.0) );
        f1s.push_back( Vec3(((x1*-1) + x2 + x3 - x4 - x5 + x6 + x7 - x8) / 8.0) );
        f2s.push_back( Vec3(((x1*-1) - x2 + x3 + x4 - x5 - x6 + x7 + x8) / 8.0) );
        f3s.push_back( Vec3(((x1*-1) - x2 - x3 - x4 + x5 + x6 + x7 + x8) / 8.0) );
        f4s.push_back( Vec3((x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8) / 8.0) );
        f5s.push_back( Vec3((x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8) / 8.0) );
        f6s.push_back( Vec3((x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8) / 8.0) );
        f7s.push_back( Vec3(((x1*-1) + x2 - x3 + x4 + x5 - x6 + x7 - x8) / 8.0) );
    }

    // store the mapping between femocs and deal.ii hexahedra
    int deal_hex_index = 0;
    map_femocs2deal = vector<int>(n_elems, -1);
    for (int i = 0; i < n_elems; ++i) {
        if (hexs->get_marker(i) > 0)
            map_femocs2deal[i] = deal_hex_index++;
    }

    map_deal2femocs = vector<int>(deal_hex_index);
    deal_hex_index = 0;
    for (int i = 0; i < n_elems; ++i) {
        if (map_femocs2deal[i] >= 0)
            map_deal2femocs[deal_hex_index++] = i;
    }
}

/* The inspiration for mapping the point was taken from
 * https://www.grc.nasa.gov/www/winddocs/utilities/b4wind_guide/trilinear.html
 */
void LinearHexahedra::project_to_nat_coords(double &u, double &v, double &w,
        const Vec3& point, const int hex) const {

    double du, dv, dw, D;
    Vec3 f0 = point - f0s[hex];
    Vec3 f1 = f1s[hex];
    Vec3 f2 = f2s[hex];
    Vec3 f3 = f3s[hex];
    Vec3 f4 = f4s[hex];
    Vec3 f5 = f5s[hex];
    Vec3 f6 = f6s[hex];
    Vec3 f7 = f7s[hex];

    // In 3D, direct calculation of uvw is very expensive (not to say impossible),
    // because system of 3 nonlinear equations should be solved.
    // More efficient is to perform Newton iterations to calculate uvw approximately.

    // loop until the desired accuracy is met or # max iterations is done
    u = 0; v = 0; w = 0;
    for (int i = 0; i < n_newton_iterations; ++i) {
        Vec3 f = (f0 - f1*u - f2*v - f3*w - f4*(u*v) - f5*(u*w) - f6*(v*w) - f7*(u*v*w));
        Vec3 fu = f1 + f4*v + f5*w + f7*(v*w);
        Vec3 fv = f2 + f4*u + f6*w + f7*(u*w);
        Vec3 fw = f3 + f5*u + f6*v + f7*(u*v);

        // solve the set of following linear equations for du, dv and dw
        //   | fu.x * du + fv.x * dv + fw.x * dw = f.x;
        //   | fu.y * du + fv.y * dv + fw.y * dw = f.y;
        //   | fu.z * du + fv.z * dv + fw.z * dw = f.z;

        D = determinant(fu, fv, fw);
        require(D != 0, "Invalid determinant: " + d2s(D));

        D  = 1.0 / D;
        du = determinant(f, fv, fw) * D;
        dv = determinant(fu, f, fw) * D;
        dw = determinant(fu, fv, f) * D;

        u += du;
        v += dv;
        w += dw;

        if (du * du + dv * dv + dw * dw < zero)
            return;
    }
}

void LinearHexahedra::project_to_nat_coords(double &u, double &v, double &w,
        int &n1, int &n2, int &n3, const int node) const {

    switch (node) {
        case 0: u=-1; v=-1; w=-1; n1=1; n2=3; n3=4; return;
        case 1: u=1;  v=-1; w=-1; n1=0; n2=2; n3=5; return;
        case 2: u=1;  v=1;  w=-1; n1=3; n2=1; n3=6; return;
        case 3: u=-1; v=1;  w=-1; n1=2; n2=0; n3=7; return;
        case 4: u=-1; v=-1; w=1;  n1=5; n2=7; n3=0; return;
        case 5: u=1;  v=-1; w=1;  n1=4; n2=6; n3=1; return;
        case 6: u=1;  v=1;  w=1;  n1=7; n2=5; n3=2; return;
        case 7: u=-1; v=1;  w=1;  n1=6; n2=4; n3=3; return;
    }
}

array<double,8> LinearHexahedra::shape_functions(const Vec3& point, const int hex) const {
    require(hex >= 0 && hex < size(), "Index out of bounds: " + d2s(hex));

    // Before calculating the shape function,
    // map the point from Cartesian xyz-coordinates to natural uvw-coordinates.
    double u, v, w;
    project_to_nat_coords(u, v, w, point, hex);

    // use natural coordinates to calculate shape functions
    return {
            (1 - u) * (1 - v) * (1 - w) / 8.0,
            (1 + u) * (1 - v) * (1 - w) / 8.0,
            (1 + u) * (1 + v) * (1 - w) / 8.0,
            (1 - u) * (1 + v) * (1 - w) / 8.0,
            (1 - u) * (1 - v) * (1 + w) / 8.0,
            (1 + u) * (1 - v) * (1 + w) / 8.0,
            (1 + u) * (1 + v) * (1 + w) / 8.0,
            (1 - u) * (1 + v) * (1 + w) / 8.0
    };
}

array<double,8> LinearHexahedra::shape_funs_dealii(const Vec3& point, const int hex) const {
    array<double,8> sf = shape_functions(point, hex);
    return {sf[0], sf[1], sf[4], sf[5], sf[3], sf[2], sf[7], sf[6]};
}

/* For more information, see the lecture materials in
 * https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
 */
array<Vec3,8> LinearHexahedra::shape_fun_grads(const Vec3& point, const int hex) const {
    require(hex >= 0 && hex < size(), "Index out of bounds: " + d2s(hex));
    double u, v, w;
    project_to_nat_coords(u, v, w, point, hex);

    // Transfer the coordinates of hex nodes into a matrix
    SimpleHex shex = get_cell(hex);
    array<Vec3,8> xyz;
    for (int i = 0; i < 8; ++i)
        xyz[i] = mesh->nodes[shex[i]];

    // calculate gradient of shape functions in uvw-space
    array<Vec3,8> dN = {
            Vec3(-(1-v)*(1-w), -(1-u)*(1-w), -(1-u)*(1-v)),
            Vec3( (1-v)*(1-w), -(1+u)*(1-w), -(1+u)*(1-v)),
            Vec3( (1+v)*(1-w),  (1+u)*(1-w), -(1+u)*(1+v)),
            Vec3(-(1+v)*(1-w),  (1-u)*(1-w), -(1-u)*(1+v)),
            Vec3(-(1-v)*(1+w), -(1-u)*(1+w),  (1-u)*(1-v)),
            Vec3( (1-v)*(1+w), -(1+u)*(1+w),  (1+u)*(1-v)),
            Vec3( (1+v)*(1+w),  (1+u)*(1+w),  (1+u)*(1+v)),
            Vec3(-(1+v)*(1+w),  (1-u)*(1+w),  (1-u)*(1+v))
    };
    for (int i = 0; i < 8; ++i)
        dN[i] *= 0.125;

    // calculate Jacobian
    array<Vec3,3> J;
    for (int k = 0; k < 8; ++k)
        for (int i = 0; i < 3; ++i)
            J[i] += xyz[k] * dN[k][i];

    // calculate determinant (and its inverse) of Jacobian
    double Jdet = determinant(J[0], J[1], J[2]);
    require(fabs(Jdet) > 1e-15, "Singular Jacobian can't be handled!");
    Jdet = 1.0 / Jdet;

    // calculate inverse of Jacobian
    array<Vec3, 3> Jinv = {
        Vec3(J[1][1]*J[2][2]-J[2][1]*J[1][2], J[2][0]*J[1][2]-J[1][0]*J[2][2], J[1][0]*J[2][1]-J[2][0]*J[1][1]),
        Vec3(J[2][1]*J[0][2]-J[0][1]*J[2][2], J[0][0]*J[2][2]-J[2][0]*J[0][2], J[2][0]*J[0][1]-J[0][0]*J[2][1]),
        Vec3(J[0][1]*J[1][2]-J[1][1]*J[0][2], J[1][0]*J[0][2]-J[0][0]*J[1][2], J[0][0]*J[1][1]-J[1][0]*J[0][1])
    };
    for (int i = 0; i < 3; ++i)
        Jinv[i] *= Jdet;

    array<Vec3,8> sfg;
    // calculate gradient of shape functions in xyz-space
    for (int k = 0; k < 8; ++k) {
//        sfg[k] = Vec3(0);
        for (int i = 0; i < 3; ++i)
            sfg[k] += Jinv[i] * dN[k][i];
    }
    return sfg;
}

array<Vec3, 8> LinearHexahedra::shape_fun_grads_dealii(const Vec3& point, const int hex) const {
    array<Vec3, 8> sfg = shape_fun_grads(point, hex);
    return {sfg[0], sfg[1], sfg[4], sfg[5], sfg[3], sfg[2], sfg[7], sfg[6]};
}

/* The gradient of shape functions in uvw space consists in case of
 * nodal data mostly of zeros. For example, for the 2nd node,
 * the uvw coordinates and the gradient look like
 *     u=1, v=-1, w=-1
 *     2*dN=[[-1     1     0     0     0     0     0     0]
 *           [ 0    -1     1     0     0     0     0     0]
 *           [ 0    -1     0     0     0     1     0     0]]
 *
 * Another important porperty of nodal data is, that only the neighbouring
 * nodes are included in Jacobian. Continuing with the example of 2nd node,
 * the neighbouring nodes and Jacobian look like
 *     n0=0, n1=2, n2=5
 *     2*J=[[ x1 - x0, y1 - y0, z1 - z0]
 *          [ x2 - x1, y2 - y1, z2 - z1]
 *          [ x5 - x1, y5 - y1, z5 - z1]]
 *
 * Those properties can be used to reduce the amount of necessary flops.
 *
 * The Jacobian can be expressed as
 *     J[i,:] = (xyz[node] - xyz[n_i]) * uvw[i]/2 .
 *
 * The explicit formula to calculate gradient of shape functions in uvw-space
 *     double dN[3][8] = {0};
 *     for (int i = 0; i < 3; ++i) {
 *         dN[i][node] = 0.5 * uvw[i];
 *         dN[i][n_i] = -0.5 * uvw[i];
 *     }
 *
 * And finally the gradient of shape functions in xyz-space
 *     sfg = array<Vec3,8>();
 *     for (int i = 0; i < 3; ++i) {
 *         sfg[node][i] = Jinv[i][0] * dN[0][node] + Jinv[i][1] * dN[1][node] + Jinv[i][2] * dN[2][node];
 *         for (int j = 0; j < 3; ++j)
 *             sfg[n_j][i] = Jinv[i][j] * dN[j][n_j];
 *     }
 *
 * As the explicit knowledge about dN is not needed,
 * last two expressions can be optimized to a single one.
 */
array<Vec3, 8> LinearHexahedra::shape_fun_grads(const int hex, const int node) const {
    require(hex >= 0 && hex < size(), "Invalid hex: " + d2s(hex));
    require(node >= 0 && node < 8, "Invalid node: " + d2s(node));

    double u, v, w; // natural coordinates
    int n1, n2, n3; // indices of nearest neighbouring nodes of the node
    project_to_nat_coords(u, v, w, n1, n2, n3, node);

    Vec3 uvw(u, v, w);
    int nnn[3] = {n1, n2, n3};
    SimpleHex shex = get_cell(hex);

    // calculate Jacobian
    Vec3 vec0 = nodes->get_vertex(shex[node]);
    array<Vec3, 3> J = {
        Vec3(vec0 - nodes->get_vertex(shex[nnn[0]])) * (uvw[0] * 0.5),
        Vec3(vec0 - nodes->get_vertex(shex[nnn[1]])) * (uvw[1] * 0.5),
        Vec3(vec0 - nodes->get_vertex(shex[nnn[2]])) * (uvw[2] * 0.5)
    };

    // calculate determinant (and its inverse) of Jacobian
    double Jdet = determinant(J[0], J[1], J[2]);
    require(fabs(Jdet) > 1e-15, "Singular Jacobian can't be handled!");
    Jdet = 1.0 / Jdet;

    // calculate inverse of Jacobian
    array<Vec3,3> Jinv = {
        Vec3(J[1][1]*J[2][2]-J[2][1]*J[1][2], J[2][1]*J[0][2]-J[0][1]*J[2][2], J[0][1]*J[1][2]-J[1][1]*J[0][2]),
        Vec3(J[2][0]*J[1][2]-J[1][0]*J[2][2], J[0][0]*J[2][2]-J[2][0]*J[0][2], J[1][0]*J[0][2]-J[0][0]*J[1][2]),
        Vec3(J[1][0]*J[2][1]-J[2][0]*J[1][1], J[2][0]*J[0][1]-J[0][0]*J[2][1], J[0][0]*J[1][1]-J[1][0]*J[0][1])
    };
    for (int i = 0; i < 3; ++i)
        Jinv[i] *= Jdet;

    // calculate gradient of shape functions in xyz-space
    array<Vec3,8> sfg;
    for (int i = 0; i < 3; ++i) {
        sfg[node][i] = 0.5 * (Jinv[i].dotProduct(uvw));
        for (int j = 0; j < 3; ++j)
            sfg[nnn[j]][i] = -0.5 * Jinv[i][j] * uvw[j];
    }

    return sfg;
}

bool LinearHexahedra::point_in_cell(const Vec3 &point, const int cell) const {
    static constexpr int n_hexs_per_tet = 4;

    int tet_index = (int) cell / n_hexs_per_tet;
    array<double,4> bcc = lintet->shape_functions(point, tet_index);

    if (bcc[0] >= 0 && bcc[1] >= 0 && bcc[2] >= 0 && bcc[3] >= 0) {
        int section = cell % n_hexs_per_tet;
        switch (section) {
            case 0:
                return bcc[0] >= bcc[1] && bcc[0] >= bcc[2] && bcc[0] >= bcc[3];
            case 1:
                return bcc[1] >= bcc[0] && bcc[1] >= bcc[2] && bcc[1] >= bcc[3];
            case 2:
                return bcc[2] >= bcc[0] && bcc[2] >= bcc[1] && bcc[2] >= bcc[3];
            case 3:
                return bcc[3] >= bcc[0] && bcc[3] >= bcc[1] && bcc[3] >= bcc[2];
        }
    }

    return false;
}

/*
 * Function uses the fact that hexahedra are always uniquely tied with the nodes of tetrahedra.
 * Therefore, knowing the barycentric coordinates inside a tetrahedron
 * makes it possible to use them to determine, in which part of the tetrahedron and therefore
 * also inside which hexahedron the point is located.
 */
int LinearHexahedra::locate_cell(const Point3 &point, const int cell_guess) const {
    static constexpr int n_hexs_per_tet = 4;

    int tet = (int) cell_guess / n_hexs_per_tet;
    tet = lintet->locate_cell(point, tet);
    int sign = 1;
    if (tet < 0) sign = -1;
    tet = abs(tet);

    // calculate barycentric coordinates for a point
    array<double,4> bcc = lintet->shape_functions(point, tet);

    // point inside a hex connected to 1st tetrahedral node ?
    if (bcc[0] >= bcc[1] && bcc[0] >= bcc[2] && bcc[0] >= bcc[3])
        return sign * (n_hexs_per_tet * tet + 0);
    // point inside a hex connected to 2nd tetrahedral node ?
    if (bcc[1] >= bcc[0] && bcc[1] >= bcc[2] && bcc[1] >= bcc[3])
        return sign * (n_hexs_per_tet * tet + 1);
    // point inside a hex connected to 3rd tetrahedral node ?
    if (bcc[2] >= bcc[0] && bcc[2] >= bcc[1] && bcc[2] >= bcc[3])
        return sign * (n_hexs_per_tet * tet + 2);
    // point inside a hex connected to 4th tetrahedral node ?
    if (bcc[3] >= bcc[0] && bcc[3] >= bcc[1] && bcc[3] >= bcc[2])
        return sign * (n_hexs_per_tet * tet + 3);

    return -1;
}

/* ==================================================================
 *  ======================= LinearTriangles ========================
 * ================================================================== */

LinearTriangles::LinearTriangles() :
        InterpolatorCells<3>(), tris(NULL), lintet(NULL) {}

LinearTriangles::LinearTriangles(const InterpolatorNodes* n, const LinearTetrahedra* lt) :
        InterpolatorCells<3>(n), tris(NULL), lintet(lt) {}

void LinearTriangles::reserve(const int n) {
    InterpolatorCells<3>::reserve(n);

    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear();  pvec.reserve(n);
    norms.clear(); norms.reserve(n);
    max_distance.clear(); max_distance.reserve(n);
}

void LinearTriangles::precompute() {
    require(mesh && tris && lintet, "NULL pointers can't be used!");
    const int n_faces = tris->size();
    const int n_nodes = mesh->nodes.size();

    expect(n_nodes > 0 && n_faces > 0, "Interpolator expects non-empty mesh!");

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Store the constant for smoothing
    decay_factor = -1.0 / tris->stat.edgemax;

    // Calculate which triangles are connected to which node
    vector<vector<int>> node2tris(n_nodes);
    for (int tri = 0; tri < n_faces; ++tri)
        for (int node : (*tris)[tri]) {
            node2tris[node].push_back(tri);
    }

    // Loop through all the faces
    for (int tri = 0; tri < n_faces; ++tri) {
        SimpleFace sface = (*tris)[tri];

        // store next nearest neighbours of a triangle
        for (int node : sface)
            for (int nbor_tri : node2tris[node]) {
                if (nbor_tri != tri)
                    neighbours[tri].push_back(nbor_tri);
            }

        Vec3 v0 = mesh->nodes.get_vec(sface[0]);
        Vec3 v1 = mesh->nodes.get_vec(sface[1]);
        Vec3 v2 = mesh->nodes.get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = tris->get_norm(tri).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);

        // calculate norms of triangles
        norms.push_back(tris->get_norm(tri));
        // store max distance from given triangle
        max_distance.push_back(e2.norm());
        // calculate centroids of triangles
        centroids.push_back(tris->get_centroid(tri));
    }
}

bool LinearTriangles::point_in_cell(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < -zero || u > 1 + zero) return false;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = qvec.dotProduct(norms[face]);
    if (v < -zero || u + v > 1 + zero) return false; // Check second & third barycentric coordinate

    // finally check the distance of the point from the triangle
    return fabs(qvec.dotProduct(edge2[face])) < max_distance[face];
}

array<double,3> LinearTriangles::shape_functions(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    const double v = tvec.dotProduct(pvec[face]);
    const double w = qvec.dotProduct(norms[face]);
    const double u = 1.0 - v - w;

    return {zero + u, zero + v, zero + w};
}

Solution LinearTriangles::interp_solution(const Point3 &point, const int t) const {
    require(mesh->tets.size() == lintet->size(),
            "Mismatch between tetrahedral mesh and interpolator sizes (" + d2s(mesh->tets.size()) + " vs " + d2s(lintet->size())
            + ")\nindicates, that LinearTetrahedra are not properly pre-computed!");

    int tri = abs(t);
    array<int, 2> tets = tris->to_tets(tri);

    require(tets[0] >= 0 && tets[1] >= 0, "Triangle " + d2s(tri) + " should have two associated tetrahedra!");

    // if the point is exactly inside the face, 3D cell locator doesn't work
    // therefore it's necessary to assume that for points on the face plane
    // either of the cells are fine
    double distance_to_tri = abs(fast_distance(point, tri));
    if (distance_to_tri <= 100.0 * zero)
        return lintet->interp_solution(point, tets[0]);

    // test if point is inside the cell directly connected to face
    if (lintet->point_in_cell(point, tets[0]))
        return lintet->interp_solution(point, tets[0]);

    // no, use tetrahedral tools to obtain the result
    int tet = lintet->locate_cell(point, tets[1]);
    return lintet->interp_solution(point, tet);
}

void LinearTriangles::interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) const {
    const int n_atoms = atoms.size();
    const int n_nodes = nodes->size();

    vector<double> bcc_sum(n_nodes); // sum of barycentric coordinates from given node
    vector<int> atom2cell(n_atoms);  // map storing the face indices that correspond to atom sequence
    scalars = vector<double>(n_atoms);
    array<double,3> weights;

    // calculate the sum of all the weights in all the nodes
    // it is neccesary to ensure that the total sum of a interpolated scalar does not change
    int cell = 0;
    for (int i = 0; i < n_atoms; ++i) {
        Point3 point = atoms[i].point;
        // Find the cell that matches best to the point
        cell = abs(this->locate_cell(point, cell));
        // calculate barycentric coordinates
        weights = shape_functions(point, cell);
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
        int tri = atom2cell[i];
        // calculate barycentric coordinates
        weights = shape_functions(atoms[i].point, tri);
        // perform interpolation
        int j = 0;
        for (int node : get_cell(tri))
            scalars[i] += nodes->get_scalar(node) * weights[j++] / bcc_sum[node];
    }
}

int LinearTriangles::near_surface(const Vec3& point, const double r_cut) const {
    require(r_cut > 0, "Invalid distance from surface: " + d2s(r_cut));

    for (int face = 0; face < size(); ++face) {
        const double dist = distance(point, face);
        if (dist >= -0.3*r_cut && dist <= r_cut) return face;
    }

    return -1;
}

double LinearTriangles::fast_distance(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    return edge2[face].dotProduct(qvec);
}

double LinearTriangles::distance(const Vec3& point, const int face) const {
    // Constants to specify the tolerances in searching outside the triangle
    const double zero = -0.1;
    const double one = 1.1;

    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < zero || u > one) return 1e100;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = norms[face].dotProduct(qvec);
    if (v < zero || u + v > one) return 1e100; // Check second & third barycentric coordinate

    // return the distance from point to triangle
    return edge2[face].dotProduct(qvec);
}

void LinearTriangles::write_vtk_cell_data(ofstream& out) const {
    InterpolatorCells<3>::write_vtk_cell_data(out);

    // write face norms
    out << "VECTORS norm double\n";
    for (int i = 0; i < size(); ++i)
        out << norms[i] << "\n";
}

/* ==================================================================
 *  ====================== QuadraticTriangles ======================
 * ================================================================== */

QuadraticTriangles::QuadraticTriangles() :
        InterpolatorCells<6>(), tris(NULL), lintri(NULL), quadtet(NULL) {}

QuadraticTriangles::QuadraticTriangles(const InterpolatorNodes* n, const LinearTriangles* lt, const QuadraticTetrahedra* qt) :
    InterpolatorCells<6>(n), tris(NULL), lintri(lt), quadtet(qt) {}

void QuadraticTriangles::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));
    markers = vector<int>(N);
    cells.clear();
    cells.reserve(N);
}

void QuadraticTriangles::precompute() {
    require(mesh && tris && lintri && quadtet, "NULL pointers can't be used!");
    const int n_faces = tris->size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Store the constant for smoothing
    decay_factor = -1.0 / tris->stat.edgemax;

    for (int i = 0; i < n_faces; ++i) {
        // calc & store triangles
        cells.push_back(calc_cell(i));
    }
}

bool QuadraticTriangles::point_in_cell(const Vec3& point, const int face) const {
    return lintri->point_in_cell(point, face);
}

int QuadraticTriangles::locate_cell(const Point3 &point, const int cell_guess) const {
    return lintri->locate_cell(point, cell_guess);
}

array<double,6> QuadraticTriangles::shape_functions(const Vec3& point, const int face) const {
    array<double,3> bcc = lintri->shape_functions(point, face);

    const double u = bcc[0];
    const double v = bcc[1];
    const double w = bcc[2];

    return {
        u * (2 * u - 1),
        v * (2 * v - 1),
        w * (2 * w - 1),
        4 * u * v,
        4 * v * w,
        4 * w * u
    };
//    return {u, v, w, 0, 0, 0};
}

Solution QuadraticTriangles::interp_solution(const Point3 &point, const int t) const {
    require(tris->size() == lintri->size(),
            "Mismatch between triangular mesh and interpolator sizes (" + d2s(tris->size()) + " vs " + d2s(lintri->size())
            + ")\nindicates, that LinearTriangles are not properly pre-computed!");

    int tri = abs(t);
    array<int, 2> tets = tris->to_tets(tri);

    require(tets[0] >= 0 && tets[1] >= 0, "Triangle " + d2s(tri) + " should have two associated tetrahedra!");

    // if the point is exactly inside the face, 3D cell locator doesn't work
    // therefore it's necessary to assume that for points on the face plane
    // either of the cells are fine
    double distance_to_tri = abs(lintri->fast_distance(point, tri));
    if (distance_to_tri <= 100.0 * zero)
        return quadtet->interp_solution(point, tets[0]);

    // test if point is inside the cell directly connected to face
    if (quadtet->point_in_cell(point, tets[0]))
        return quadtet->interp_solution(point, tets[0]);

    // no, use tetrahedral tools to obtain the result
    int tet = quadtet->locate_cell(point, tets[1]);
    return quadtet->interp_solution(point, tet);
}

SimpleCell<6> QuadraticTriangles::calc_cell(const int tri) const {
    require(tri >= 0 && tri < tris->size(), "Invalid index: " + d2s(tri));
    if (mesh->quads.size() == 0)
        return QuadraticTri(0);

    static constexpr int n_quads_per_tri = 3;
    array<vector<unsigned>, n_quads_per_tri> edge_nodes;

    // locate hexahedral nodes that are located in the middle of edges
    for (int i = 0; i < n_quads_per_tri; ++i) {
        for (int quadnode : mesh->quads[n_quads_per_tri * tri + i])
            if (mesh->nodes.get_marker(quadnode) == TYPES.EDGECENTROID)
                edge_nodes[i].push_back(quadnode);
    }

    // find second order nodes
    const int n4 = common_entry(edge_nodes[0], edge_nodes[1]);
    const int n5 = common_entry(edge_nodes[1], edge_nodes[2]);
    const int n6 = common_entry(edge_nodes[2], edge_nodes[0]);

    return QuadraticTri((*tris)[tri], n4, n5, n6);
}

/* ==================================================================
 *  ====================== LinearQuadrangles =======================
 * ================================================================== */

LinearQuadrangles::LinearQuadrangles() :
        InterpolatorCells<4>(), quads(NULL), lintri(NULL), linhex(NULL) {}

LinearQuadrangles::LinearQuadrangles(const InterpolatorNodes* n, const LinearTriangles* lt, const LinearHexahedra* lh) :
        InterpolatorCells<4>(n), quads(NULL), lintri(lt), linhex(lh) {}

void LinearQuadrangles::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + d2s(N));
    markers = vector<int>(N);
}

void LinearQuadrangles::precompute() {
    require(mesh && quads && lintri && linhex, "NULL pointers can't be used!");
    require(mesh->tris.size() == lintri->size(),
            "Mismatch between triangular mesh and interpolator sizes (" + d2s(mesh->tris.size()) + " vs " + d2s(lintri->size())
            + ")\nindicates, that LinearTriangles are not properly pre-computed!");


    const int n_faces = quads->size();
    expect(n_faces > 0, "Interpolator expects non-empty mesh!");
    reserve(n_faces);

//    for (int i = 0; i < n_faces; ++i) {
//        // make the markers to correspond to lintri
//        markers[i] = lintri->get_marker(quads->to_tri(i));
//    }
}

Solution LinearQuadrangles::interp_solution(const Point3 &point, const int q) const {
    int quad = abs(q);
    array<int, 2> hexs = quads->to_hexs(quad);

    require(hexs[0] >= 0 && hexs[1] >= 0, "Quadrangle " + d2s(quad) + " should have two associated hexahedra!");

    // if the point is exactly inside the face, 3D cell locator doesn't work
    double distance_to_tri = abs( lintri->fast_distance(point, mesh->quads.to_tri(quad)) );
    if (distance_to_tri <= 100.0 * zero)
        return linhex->interp_solution(point, hexs[0]);

    // test if point is inside the 3D cell that is associated with the face
    if (linhex->point_in_cell(point, hexs[0]))
        return linhex->interp_solution(point, hexs[0]);

    // no, use hexahedral tools to obtain the result
    int hex = linhex->locate_cell(point, hexs[1]);
    return linhex->interp_solution(point, hex);
}

bool LinearQuadrangles::point_in_cell(const Vec3 &point, const int cell) const {
    static constexpr int n_quads_per_tri = 3;

    int tri = (int) cell / n_quads_per_tri;
    array<double,3> bcc = lintri->shape_functions(point, tri);

    if (bcc[0] >= 0 && bcc[1] >= 0 && bcc[2] >= 0) {
        int section = cell % n_quads_per_tri;
        switch (section) {
            case 0:
                return bcc[0] >= bcc[1] && bcc[0] >= bcc[2];
            case 1:
                return bcc[1] >= bcc[0] && bcc[1] >= bcc[2];
            case 2:
                return bcc[2] >= bcc[0] && bcc[2] >= bcc[1];
        }
    }

    return false;
}

/*
 * Function uses the fact that hexahedra are always uniquely tied with the nodes of tetrahedra.
 * Therefore, knowing the barycentric coordinates inside a tetrahedron
 * makes it possible to use them to determine, in which part of the tetrahedron and therefore
 * also inside which hexahedron the point is located.
 */
int LinearQuadrangles::locate_cell(const Point3 &point, const int cell_guess) const {
    int tri = quads->to_tri(cell_guess);
    tri = lintri->locate_cell(point, tri);
    int sign = 1;
    if (tri < 0) sign = -1;
    tri = abs(tri);

    // list of quadrangles where the point might be located
    array<int,3> quad_indices = mesh->tris.to_quads(tri);

    // calculate barycentric coordinates for a point
    array<double,3> bcc = lintri->shape_functions(point, tri);

    // point inside a quad connected to 1st triangular node ?
    if (bcc[0] >= bcc[1] && bcc[0] >= bcc[2])
        return sign * quad_indices[0];
    // point inside a quad connected to 2nd triangular node ?
    if (bcc[1] >= bcc[0] && bcc[1] >= bcc[2])
        return sign * quad_indices[1];
    // point inside a quad connected to 3rd triangular node ?
    if (bcc[2] >= bcc[0] && bcc[2] >= bcc[1])
        return sign * quad_indices[2];

    return -1;
}

template class InterpolatorCells<3> ;
template class InterpolatorCells<6> ;
template class InterpolatorCells<4> ;
template class InterpolatorCells<10>;
template class InterpolatorCells<8>;

} /* namespace femocs */
