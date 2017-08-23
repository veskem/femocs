/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "Macros.h"
#include "Media.h"
#include "LinearInterpolator.h"

#include <float.h>

using namespace std;
namespace femocs {

/* ==================================================================
 *  ===================== LinearInterpolator =======================
 * ================================================================== */

// Initialize data vectors
LinearInterpolator::LinearInterpolator() : mesh(NULL) {
    reserve(0);
    reserve_precompute(0);
}

LinearInterpolator::LinearInterpolator(const TetgenMesh* m) : mesh(m) {
    reserve(0);
    reserve_precompute(0);
}

// Add solution to solutions vector
void LinearInterpolator::append_solution(const Solution& solution) {
    expect((unsigned)size() < solutions.capacity(), "Allocated vector size exceeded!");
    solutions.push_back(solution);
}

// Enable or disable the search of points slightly outside the tetrahedron
void LinearInterpolator::search_outside(const bool enable) {
    if (enable) {
        zero = -1.0 * epsilon;
        one = 1.0 + epsilon;
    } else {
        zero = 0;
        one = 1.0;
    }
}

// Reserve memory for interpolation data
void LinearInterpolator::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + to_string(N));
    atoms.clear();
    solutions.clear();

    atoms.reserve(N);
    solutions.reserve(N);
}

// Reserve memory for pre-compute data
void LinearInterpolator::reserve_precompute(const int N) {
    neighbours = vector<vector<int>>(N);
    centroids.reserve(N);
}

// Return full solution on i-th node
Solution LinearInterpolator::get_solution(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solutions[i];
}

// Return vector component of solution on i-th node
Vec3 LinearInterpolator::get_vector(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solutions[i].vector;
}

// Return scalar component of solution on i-th node
double LinearInterpolator::get_scalar(const int i) const {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solutions[i].scalar;
}

// Compile data string from the data vectors for file output
string LinearInterpolator::get_data_string(const int i) const {
    if (i < 0)
        return "LinearInterpolator properties=id:I:1:pos:R:3:marker:I:1:force:R:3:elfield_norm:R:1:potential:R:1";

    ostringstream strs;
    strs << atoms[i] << " " << solutions[i];
    return strs.str();
}

void LinearInterpolator::get_cell_data(ofstream& out) const {
    const int celltype = get_cell_type();
    const int n_cells = get_n_cells();
    const int n_atoms = size();

    get_cells(out);

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";

    out << "\nPOINT_DATA " << n_atoms << "\n";

    // write atom IDs
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].id << "\n";

    // write atom markers
    out << "SCALARS marker int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << atoms[i].marker << "\n";

    out << "SCALARS elfield_norm double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << solutions[i].norm << "\n";

    out << "SCALARS potential double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n_atoms; ++i)
        out << solutions[i].scalar << "\n";
}

// Determinant of 3x3 matrix which's last column consists of ones
double LinearInterpolator::determinant(const Vec3 &v1, const Vec3 &v2) {
    return v1.x * (v2.y - v2.z) - v1.y * (v2.x - v2.z) + v1.z * (v2.x - v2.y);
}

// Determinant of 3x3 matrix which's columns consist of Vec3-s
double LinearInterpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) {
    return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z)
            + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

// Determinant of 4x4 matrix which's last column consists of ones
double LinearInterpolator::determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3,
        const Vec3 &v4) {
    const double det1 = determinant(v2, v3, v4);
    const double det2 = determinant(v1, v3, v4);
    const double det3 = determinant(v1, v2, v4);
    const double det4 = determinant(v1, v2, v3);

    return det4 - det3 + det2 - det1;
}

// Determinant of 4x4 matrix which's columns consist of Vec4-s
double LinearInterpolator::determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3,
        const Vec4 &v4) {
    double det1 = determinant(Vec3(v1.y, v1.z, v1.w), Vec3(v2.y, v2.z, v2.w),
            Vec3(v3.y, v3.z, v3.w));
    double det2 = determinant(Vec3(v1.x, v1.z, v1.w), Vec3(v2.x, v2.z, v2.w),
            Vec3(v3.x, v3.z, v3.w));
    double det3 = determinant(Vec3(v1.x, v1.y, v1.w), Vec3(v2.x, v2.y, v2.w),
            Vec3(v3.x, v3.y, v3.w));
    double det4 = determinant(Vec3(v1.x, v1.y, v1.z), Vec3(v2.x, v2.y, v2.z),
            Vec3(v3.x, v3.y, v3.z));

    return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
}

/* ==================================================================
 *  =================== TetrahedronInterpolator ====================
 * ================================================================== */

// Initialize data vectors
TetrahedronInterpolator::TetrahedronInterpolator() : LinearInterpolator(),
        radius1(0), radius2(0), E0(0) {}

TetrahedronInterpolator::TetrahedronInterpolator(const TetgenMesh* m) : LinearInterpolator(m),
        radius1(0), radius2(0), E0(0) {}

// Set parameters for calculating analytical solution
void TetrahedronInterpolator::set_analyt(const Point3& origin, const double E0, const double radius1, const double radius2) {
    this->E0 = E0;
    this->radius1 = radius1;
    if (radius2 > radius1)
        this->radius2 = radius2;
    else
        this->radius2 = radius1;
    this->origin = origin;
}

// Analytical potential for i-th point near the hemisphere
double TetrahedronInterpolator::get_analyt_potential(const int i, const Point3& origin) const {
    Point3 point = get_atom(i).point - origin;
    double r = point.distance(Point3(0));
    return -E0 * point.z * (1 - pow(radius1 / r, 3.0));
}

// Analytical electric field for i-th point near the hemisphere
Vec3 TetrahedronInterpolator::get_analyt_field(const int i) const {
    Point3 point = get_atom(i).point - origin;
    double r5 = pow(point.x * point.x + point.y * point.y + point.z * point.z, 2.5);
    double r3 = pow(radius1, 3.0);
    double f = point.x * point.x + point.y * point.y - 2.0 * point.z * point.z;

    double Ex = 3 * E0 * r3 * point.x * point.z / r5;
    double Ey = 3 * E0 * r3 * point.y * point.z / r5;
    double Ez = E0 * (1.0 - r3 * f / r5);

    return Vec3(Ex, Ey, Ez);
}

double TetrahedronInterpolator::get_enhancement() const {
    double Emax = -error_field;
    for (Solution s : solutions)
        if (s.norm != error_field && s.norm > Emax) Emax = s.norm;

    return fabs(Emax / E0);
}

double TetrahedronInterpolator::get_analyt_enhancement() const {
    expect(radius1 > 0, "Invalid minor semi-axis: " + to_string(radius1));

    if ( radius2 <= radius1 )
        return 3.0;
    else {
        double nu = radius2 / radius1;
        double zeta = sqrt(nu*nu - 1);
        return pow(zeta, 3.0) / (nu * log(zeta + nu) - zeta);
    }
}

void TetrahedronInterpolator::print_enhancement() const {
    double gamma1 = get_enhancement();
    double gamma2 = get_analyt_enhancement();

    stringstream stream;
    stream << fixed << setprecision(3);
    stream << "field enhancements:  Femocs:" << gamma1
            << "  analyt:" << gamma2
            << "  f-a:" << gamma1 - gamma2
            << "  f/a:" << gamma1 / gamma2;

    write_verbose_msg(stream.str());
}

// Print the deviation from the analytical solution of hemi-ellipsoid on the infinite surface
void TetrahedronInterpolator::print_error(const Coarseners& c) const {
    if (!MODES.VERBOSE) return;

    double rms_error = 0;
    int n_points = 0;
    for (int i = 0; i < size(); ++i)
        // look for the tetrahedral nodes whose potential == 0, i.e that are on the surface
        if ( get_marker(i) == TYPES.TETNODE && solutions[i].scalar == 0.0 &&
                c.inside_interesting_region(get_point(i)) )
        {
            double analyt = get_analyt_field(i).norm();
            double numerical = solutions[i].norm;
//            cout << get_point(i) << " " << analyt << " " << numerical << endl;

            double error = (numerical - analyt) / numerical;
            rms_error += error * error;
            n_points++;
        }
    require(n_points > 0, "No tetrahedral points on the surface of nanotip!");
    rms_error = sqrt(rms_error / n_points);

    stringstream stream;
    stream << fixed << setprecision(3);
    stream << "rms surface error: " << rms_error;

    write_verbose_msg(stream.str());
}

// Print statistics about solution on node points
void TetrahedronInterpolator::print_statistics() const {
    if (!MODES.VERBOSE) return;

    const int n_atoms = size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        double s = solutions[i].scalar;
        Vec3 v = solutions[i].vector;
        if (s >= error_field) continue;

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

/* Return the mapping between tetrahedral & hexahedral mesh nodes,
 nodes & hexahedral elements and nodes & element's vertices  */
void TetrahedronInterpolator::get_maps(vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs,
        dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh, const double eps) {

    require(tria->n_vertices() > 0, "Invalid triangulation size!");

    const int n_femocs_nodes = size();
    const int n_dealii_nodes = tria->n_used_vertices();
    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;

    vector<int> node2hex(n_dealii_nodes), node2vert(n_dealii_nodes);
    tet2hex = vector<int>(n_femocs_nodes, -1);

    typename dealii::Triangulation<3>::active_vertex_iterator vertex = tria->begin_active_vertex();
    // Loop through tetrahedral mesh vertices
    for (int i = 0; i < n_femocs_nodes && vertex != tria->end_vertex(); ++i)
        if ( get_point(i).distance2(vertex->vertex()) < eps ) {
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

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
bool TetrahedronInterpolator::average_tetnodes() {
//    return 0;

    const double decay_factor = -1.0 / mesh->elems.stat.edgemax;

    vector<vector<unsigned int>> cells;
    mesh->calc_pseudo_vorocells(cells);

    bool fail = false;

    // loop through the tetrahedral nodes
    for (int i = 0; i < mesh->nodes.stat.n_tetnode; ++i) {
        if (solutions[i].norm >= error_field) continue;

        Point3 tetnode = mesh->nodes[i];
        Vec3 vec(0);
        double w_sum = 0;

        // tetnode new solution will be the weighed average of the solutions on its voronoi cell nodes
        for (unsigned v : cells[i]) {
            double w = exp(decay_factor * tetnode.distance(mesh->nodes[v]));
            w_sum += w;
            vec += solutions[v].vector * w;
        }

        if (w_sum > 0) {
            solutions[i].vector = vec * (1.0 / w_sum);
            solutions[i].norm = solutions[i].vector.norm();
        } else {
            expect(false, "Tetrahedral node " + to_string(i) + " can't be averaged!");
            fail = true;
        }
    }

    return fail;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::Laplace<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = mesh->nodes.size();
    const double eps = 1e-5 * mesh->elems.stat.edgemin;

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int i = 0; i < n_nodes; ++i)
        append( Atom(i, mesh->nodes[i], mesh->nodes.get_marker(i)) );

    // Precompute tetrahedra to make interpolation faster
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
            append_solution(Solution(error_field));
    }

    // remove the spikes in the solution
    if (average_tetnodes())
        return true;

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : mesh->elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
                return true;

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::CurrentsAndHeatingStationary<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = mesh->nodes.size();
    const double eps = 1e-5 * mesh->elems.stat.edgemin;

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int n = 0; n < n_nodes; ++n)
        append( Atom(n, mesh->nodes[n], mesh->nodes.get_marker(n)) );

    // Precompute tetrahedra to make interpolation faster
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
            append_solution(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : mesh->elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
                return true;

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool TetrahedronInterpolator::extract_solution(fch::CurrentsAndHeating<3>* fem) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = mesh->nodes.size();
    const double eps = 1e-5 * mesh->elems.stat.edgemin;

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int n = 0; n < n_nodes; ++n)
        append( Atom(n, mesh->nodes[n], mesh->nodes.get_marker(n)) );

    // Precompute tetrahedra to make interpolation faster
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
            append_solution(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : mesh->elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
                return true;

    return false;
}

// Reserve memory for pre-compute data
void TetrahedronInterpolator::reserve_precompute(const int N) {
    LinearInterpolator::reserve_precompute(N);

    tetrahedra.clear();
    det0.clear();
    det1.clear();
    det2.clear();
    det3.clear();
    det4.clear();
    tet_not_valid.clear();

    tetrahedra.reserve(N);
    det0.reserve(N);
    det1.reserve(N);
    det2.reserve(N);
    det3.reserve(N);
    det4.reserve(N);
    tet_not_valid.reserve(N);
}

// Precompute the data to tetrahedra to make later bcc calculations faster
void TetrahedronInterpolator::precompute() {
    const int n_elems = mesh->elems.size();
    reserve_precompute(n_elems);
    double d0, d1, d2, d3, d4;

    // Copy tetrahedra
    for (int i = 0; i < n_elems; ++i)
        tetrahedra.push_back(mesh->elems[i]);

    // Copy tetrahedra neighbours
    for (int i = 0; i < n_elems; ++i)
        neighbours.push_back(mesh->elems.get_neighbours(i));

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroids.push_back(mesh->elems.get_centroid(i));

    /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
     * barycentric coordinate of tetrahedra using the relations below */
    for (SimpleElement se : tetrahedra) {
        Vec3 v1 = mesh->nodes.get_vec(se[0]);
        Vec3 v2 = mesh->nodes.get_vec(se[1]);
        Vec3 v3 = mesh->nodes.get_vec(se[2]);
        Vec3 v4 = mesh->nodes.get_vec(se[3]);

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
bool TetrahedronInterpolator::point_in_tetrahedron(const Point3 &point, const int i) {
    require(i >= 0 && i < static_cast<int>(det0.size()), "Index out of bounds: " + to_string(i));

    // Ignore co-planar tetrahedra
    // no need to check because Tetgen guarantees non-co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < zero) return false;

    // All bcc-s are >= 0, so point is inside the tetrahedron
    return true;
}

// Calculate barycentric coordinates for point
Vec4 TetrahedronInterpolator::get_bcc(const Point3 &point, const int elem) const {
    require(elem >= 0 && elem < static_cast<int>(det0.size()),
            "Index out of bounds: " + to_string(elem));

    Vec4 pt(point, 1);
    double bcc1 = det0[elem] * pt.dotProduct(det1[elem]);
    double bcc2 = det0[elem] * pt.dotProduct(det2[elem]);
    double bcc3 = det0[elem] * pt.dotProduct(det3[elem]);
    double bcc4 = det0[elem] * pt.dotProduct(det4[elem]);

    return Vec4(bcc1, bcc2, bcc3, bcc4);
}

// Find the element which contains the point or is the closest to it
int TetrahedronInterpolator::locate_element(const Point3 &point, const int elem_guess) {
    // Check the guessed element
    if (point_in_tetrahedron(point, elem_guess)) return elem_guess;

    const int n_elems = det0.size();
    const int n_nbor_layers = 6;  // choose the amount of nearest neighboring layers that are checked before the full search

    vector<vector<int>> nbors(n_nbor_layers);
    vector<bool> elem_checked(n_elems);

    elem_checked[elem_guess] = true;

    // Check all tetrahedra on the given neighbouring layer
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of neighbour list
        if (layer == 0)
            nbors[0] = neighbours[elem_guess];
        else {
            for (int nbor : nbors[layer-1])
                if (nbor >= 0)
                    nbors[layer].insert(nbors[layer].end(), neighbours[nbor].begin(), neighbours[nbor].end());
        }

        // check whether some of the unchecked neighbouring tetrahedra surround the point
        for (int elem : nbors[layer])
            if (elem >= 0 && !elem_checked[elem]) {
                if (point_in_tetrahedron(point, elem))
                    return elem;
                else
                    elem_checked[elem] = true;
            }
    }

    // If no success, loop through all the elements
    double min_distance2 = DBL_MAX;
    int min_index = 0;

    for (int elem = 0; elem < n_elems; ++elem) {
        // If correct element is found, we're done
        if (!elem_checked[elem] && point_in_tetrahedron(point, elem))
            return elem;

        // Otherwise look for the element whose centroid is closest to the point
        else {
            double distance2 = point.distance2(centroids[elem]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = elem;
            }
        }
    }

    // If no perfect element found, return the best.
    // Indicate the imperfectness with the minus sign
    return -min_index;
}

// Calculate interpolation or all the data for point inside or near the elem-th tetrahedron
Solution TetrahedronInterpolator::interp_solution(const Point3 &point, const int elem) const {
    require(elem >= 0 && elem < static_cast<int>(tetrahedra.size()),
            "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);
    SimpleElement selem = tetrahedra[elem];

    // Interpolate electric field
    Vec3 elfield_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        elfield_i += solutions[selem[i]].vector * bcc[i];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        potential_i += solutions[selem[i]].scalar * bcc[i];

    return Solution(elfield_i, potential_i);
}

// Calculate interpolation for vector data for point inside or near the elem-th tetrahedron
Vec3 TetrahedronInterpolator::interp_vector(const Point3 &point, const int elem) const {
    require(elem >= 0 && elem < static_cast<int>(tetrahedra.size()),
            "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);
    SimpleElement selem = tetrahedra[elem];

    // Interpolate electric field
    Vec3 vec(0.0);
    for (int i = 0; i < selem.size(); ++i)
        vec += solutions[selem[i]].vector * bcc[i];

    return vec;
}

// Calculate interpolation for scalar data for point inside or near the elem-th tetrahedron
double TetrahedronInterpolator::interp_scalar(const Point3 &point, const int elem) const {
    require(elem >= 0 && elem < static_cast<int>(tetrahedra.size()),
            "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);
    SimpleElement selem = tetrahedra[elem];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        potential_i += solutions[selem[i]].scalar * bcc[i];

    return potential_i;
}

void TetrahedronInterpolator::get_cells(ofstream& out) const {
    const int dim = tetrahedra[0].size();          // number of vertices in the cell
    const int n_cells = get_n_cells();

    // Output the vertices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << tetrahedra[i] << "\n";
}
/* ==================================================================
 *  ===================== TriangleInterpolator ======================
 * ================================================================== */

TriangleInterpolator::TriangleInterpolator() : LinearInterpolator(NULL) {}

TriangleInterpolator::TriangleInterpolator(const TetgenMesh* m) : LinearInterpolator(m) {}

// Reserve memory for precompute data
void TriangleInterpolator::reserve_precompute(const int n) {
    LinearInterpolator::reserve_precompute(n);

    edge1.clear(); edge1.reserve(n);
    edge2.clear(); edge2.reserve(n);
    vert0.clear(); vert0.reserve(n);
    pvec.clear();  pvec.reserve(n);
    is_parallel.clear(); is_parallel.reserve(n);

    triangles.clear(); triangles.reserve(n);
}

// Precompute the data needed to calculate the distance of points from surface in the direction of triangle norms
void TriangleInterpolator::precompute() {
    const int n_faces = mesh->faces.size();

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace sface = mesh->faces[i];

        Vec3 v0 = mesh->nodes.get_vec(sface[0]);
        Vec3 v1 = mesh->nodes.get_vec(sface[1]);
        Vec3 v2 = mesh->nodes.get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = mesh->faces.get_norm(i).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);
        centroids.push_back(mesh->faces.get_centroid(i));
        triangles.push_back(sface);

        // calculate the neighbour list for triangles
        for (int j = i+1; j < n_faces; ++j)
            if (sface.neighbor(mesh->faces[j])) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
    }
}

// Find the element which contains the point or is the closest to it
int TriangleInterpolator::locate_face(const Vec3 &point, const int face_guess) {
    // Check the guessed element
    if (projection_in_triangle(point, face_guess)) return face_guess;

    const int n_faces = vert0.size();
    const int n_nbor_layers = 6;  // choose the amount of nearest neighbouring layers that are checked before the full search

    vector<vector<int>> nbors(n_nbor_layers);
    vector<bool> face_checked(n_faces);

    face_checked[face_guess] = true;

    // Check all triangles on the given neighbouring layer
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of neighbour list
        if (layer == 0)
            nbors[0] = neighbours[face_guess];
        else {
            for (unsigned nbor : nbors[layer-1])
                if (nbor >= 0)
                    nbors[layer].insert(nbors[layer].end(), neighbours[nbor].begin(), neighbours[nbor].end());
        }

        // check whether some of the unchecked neighbouring triangles surround the point
        for (unsigned face : nbors[layer])
            if (face >= 0 && !face_checked[face]) {
                if (projection_in_triangle(point, face))
                    return face;
                else
                    face_checked[face] = true;
            }
    }

    // If no success, loop through all the elements
    double min_distance2 = DBL_MAX;
    int min_index = 0;
    Point3 p(point.x, point.y, point.z);

    for (int face = 0; face < n_faces; ++face) {
        // If correct face is found, we're done
        if (!face_checked[face] && projection_in_triangle(point, face))
            return face;

        // Otherwise look for the face whose centroid is closest to the point
        else {
            double distance2 = p.distance2(centroids[face]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = face;
            }
        }
    }

    // If no perfect element found, return the best.
    // Indicate the imperfectness with the minus sign
    return -min_index;
}

// Calculate barycentric coordinates for a projection of a point inside the triangle
Vec3 TriangleInterpolator::get_bcc(const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double u = tvec.dotProduct(pvec[face]);
    double v = qvec.dotProduct(mesh->faces.get_norm(face));

    return Vec3(u, v, 1.0 - u - v);
}

// Check whether the projection of a point is inside the triangle
bool TriangleInterpolator::projection_in_triangle(const Vec3& point, const int face) {
    require(face >= 0 && face < static_cast<int>(vert0.size()), "Index out of bounds: " + to_string(face));

    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < zero || u > one) return false;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = qvec.dotProduct(mesh->faces.get_norm(face));
    if (v < zero || u + v > one) return false; // Check second & third barycentric coordinate

    return true;
}

// Calculate interpolation or all the data for point inside or near the elem-th tetrahedron
Solution TriangleInterpolator::interp_solution(const Point3 &point, const int face) const {
    require(face >= 0 && face < static_cast<int>(triangles.size()),
            "Index out of bounds: " + to_string(face));

    // Get barycentric coordinates of point in tetrahedron
    Vec3 bcc = get_bcc(Vec3(point), face);
    SimpleFace cell = triangles[face];

    // Interpolate electric field
    Vec3 elfield_i(0.0);
    for (int i = 0; i < cell.size(); ++i)
        elfield_i += solutions[cell[i]].vector * bcc[i];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < cell.size(); ++i)
        potential_i += solutions[cell[i]].scalar * bcc[i];

    return Solution(elfield_i, potential_i);
}

// Calculate interpolation for vector data for point inside or near the elem-th tetrahedron
Vec3 TriangleInterpolator::interp_vector(const Point3 &point, const int face) const {
    require(face >= 0 && face < static_cast<int>(triangles.size()),
            "Index out of bounds: " + to_string(face));

    // Get barycentric coordinates of point in tetrahedron
    Vec3 bcc = get_bcc(Vec3(point), face);
    SimpleFace cell = triangles[face];

    // Interpolate electric field
    Vec3 vec(0.0);
    for (int i = 0; i < cell.size(); ++i)
        vec += solutions[cell[i]].vector * bcc[i];

    return vec;
}

// Calculate interpolation for scalar data for point inside or near the elem-th tetrahedron
double TriangleInterpolator::interp_scalar(const Point3 &point, const int elem) const {
    require(elem >= 0 && elem < static_cast<int>(triangles.size()),
            "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec3 bcc = get_bcc(Vec3(point), elem);
    SimpleFace cell = triangles[elem];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < cell.size(); ++i)
        potential_i += solutions[cell[i]].scalar * bcc[i];

    return potential_i;
}

void TriangleInterpolator::get_cells(ofstream& out) const {
    const int dim = triangles[0].size();  // number of vertices in the cell
    const int n_cells = get_n_cells();

    // Output the vertices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << triangles[i] << "\n";
}

} // namespace femocs
