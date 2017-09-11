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

// Analytical potential for i-th point near the hemisphere
template<int dim>
double LinearInterpolator<dim>::get_analyt_potential(const int i, const Point3& origin) const {
    require(i >= 0 && i < nodes->size(), "Invalid index: " + to_string(i));

    Point3 point = (*nodes)[i] - origin;
    double r = point.distance(Point3(0));
    return -E0 * point.z * (1 - pow(radius1 / r, 3.0));
}

// Analytical electric field for i-th point near the hemisphere
template<int dim>
Vec3 LinearInterpolator<dim>::get_analyt_field(const int i) const {
    require(i >= 0 && i < nodes->size(), "Invalid index: " + to_string(i));

    Point3 point = (*nodes)[i] - origin;
    double r5 = pow(point.x * point.x + point.y * point.y + point.z * point.z, 2.5);
    double r3 = pow(radius1, 3.0);
    double f = point.x * point.x + point.y * point.y - 2.0 * point.z * point.z;

    double Ex = 3 * E0 * r3 * point.x * point.z / r5;
    double Ey = 3 * E0 * r3 * point.y * point.z / r5;
    double Ez = E0 * (1.0 - r3 * f / r5);

    return Vec3(Ex, Ey, Ez);
}

// Analytical field enhancement for ellipsoidal nanotip
template<int dim>
double LinearInterpolator<dim>::get_analyt_enhancement() const {
    expect(radius1 > 0, "Invalid nanotip minor semi-axis: " + to_string(radius1));

    if ( radius2 <= radius1 )
        return 3.0;
    else {
        double nu = radius2 / radius1;
        double zeta = sqrt(nu*nu - 1);
        return pow(zeta, 3.0) / (nu * log(zeta + nu) - zeta);
    }
}

// Compare the analytical and calculated field enhancement
template<int dim>
bool LinearInterpolator<dim>::compare_enhancement() const {
    double Emax = -error_field;
    for (Solution s : solutions)
        if (s.norm != error_field && s.norm > Emax) Emax = s.norm;

    const double gamma1 = fabs(Emax / E0);
    const double gamma2 = get_analyt_enhancement();
    const double beta = fabs(gamma1 / gamma2);

    stringstream stream;
    stream << fixed << setprecision(3);
    stream << "field enhancements:  (F)emocs:" << gamma1
            << "  (A)nalyt:" << gamma2
            << "  F-A:" << gamma1 - gamma2
            << "  F/A:" << gamma1 / gamma2;

    write_verbose_msg(stream.str());
    return beta < beta_min || beta > beta_max;
}

// Interpolate both scalar and vector data inside or near the cell
template<int dim>
Solution LinearInterpolator<dim>::interp_solution(const Point3 &point, const int cell) const {
    require(cell >= 0 && cell < cells.size(), "Index out of bounds: " + to_string(cell));

    // Get barycentric coordinates of point in tetrahedron
    array<double,dim> bcc = get_bcc(Vec3(point), cell);
    SimpleCell<dim> scell = cells[cell];

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
Vec3 LinearInterpolator<dim>::interp_vector(const Point3 &point, const int cell) const {
    require(cell >= 0 && cell < cells.size(), "Index out of bounds: " + to_string(cell));

    // Get barycentric coordinates of point in tetrahedron
    array<double,dim> bcc = get_bcc(Vec3(point), cell);
    SimpleCell<dim> scell = cells[cell];

    // Interpolate electric field
    Vec3 vector_i(0.0);
    for (int i = 0; i < dim; ++i)
        vector_i += solutions[scell[i]].vector * bcc[i];

    return vector_i;
}

// Interpolate scalar data inside or near the cell
template<int dim>
double LinearInterpolator<dim>::interp_scalar(const Point3 &point, const int cell) const {
    require(cell >= 0 && cell < cells.size(), "Index out of bounds: " + to_string(cell));

    // calculate barycentric coordinates
    array<double,dim> bcc = get_bcc(Vec3(point), cell);
    SimpleCell<dim> scell = cells[cell];

    // Interpolate potential
    double scalar_i(0.0);
    for (int i = 0; i < dim; ++i)
        scalar_i += solutions[scell[i]].scalar * bcc[i];

    return scalar_i;
}

// Find the cell which contains the point or is the closest to it
template<int dim>
int LinearInterpolator<dim>::locate_cell(const Point3 &point, const int cell_guess) {
    // Check the guessed element
    Vec3 vec_point(point);
    if (point_in_cell(vec_point, cell_guess)) return cell_guess;

    const int n_cells = cells.size();
    const int n_nbor_layers = 6;  // choose the amount of nearest neighbouring layers that are checked before the full search

    vector<vector<int>> nbors(n_nbor_layers);
    vector<bool> cell_checked(n_cells);

    cell_checked[cell_guess] = true;

    // Check all triangles on the given neighbouring layer
    for (int layer = 0; layer < n_nbor_layers; ++layer) {
        // build next layer of neighbour list
        if (layer == 0)
            nbors[0] = neighbours[cell_guess];
        else {
            for (unsigned nbor : nbors[layer-1])
                if (nbor >= 0)
                    nbors[layer].insert(nbors[layer].end(), neighbours[nbor].begin(), neighbours[nbor].end());
        }

        // check whether some of the unchecked neighbouring cells surround the point
        for (unsigned cell : nbors[layer])
            if (cell >= 0 && !cell_checked[cell]) {
                if (point_in_cell(vec_point, cell))
                    return cell;
                else
                    cell_checked[cell] = true;
            }
    }

    // If no success, loop through all the elements
    double min_distance2 = 1e100;
    int min_index = 0;

    for (int cell = 0; cell < n_cells; ++cell) {
        // If correct face is found, we're done
        if (!cell_checked[cell] && point_in_cell(vec_point, cell))
            return cell;

        // Otherwise look for the face whose centroid is closest to the point
        else {
            double distance2 = point.distance2(centroids[cell]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = cell;
            }
        }
    }

    // If no perfect element found, return the best.
    // Indicate the imperfectness with the minus sign
    return -min_index;
}

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
template<int dim>
bool LinearInterpolator<dim>::average_sharp_nodes(const vector<vector<unsigned>>& vorocells,
        const double edgemax) {
    const double decay_factor = -1.0 / edgemax;

    // loop through the tetrahedral nodes
    for (int i = 0; i < vorocells.size(); ++i) {
        if (solutions[i].norm >= error_field) continue;

        Point3 tetnode = (*nodes)[i];
        Vec3 vec(0);
        double w_sum = 0;

        // tetnode new solution will be the weighed average of the solutions on its voronoi cell nodes
        for (unsigned v : vorocells[i]) {
            double w = exp(decay_factor * tetnode.distance((*nodes)[v]));
            w_sum += w;
            vec += solutions[v].vector * w;
        }

        if (w_sum > 0) {
            solutions[i].vector = vec * (1.0 / w_sum);
            solutions[i].norm = solutions[i].vector.norm();
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

    const int n_nodes = nodes->size();
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
    out << "SCALARS ID int\nLOOKUP_TABLE default\n";
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
}

/* ==================================================================
 *  =================== TetrahedronInterpolator ====================
 * ================================================================== */

// Initialize data vectors
TetrahedronInterpolator::TetrahedronInterpolator(const TetgenMesh* m) :
        LinearInterpolator<4>(m), elems(&m->elems) {}

// Print the deviation from the analytical solution of hemi-ellipsoid on the infinite surface
void TetrahedronInterpolator::print_error(const Coarseners& c) const {
    if (!MODES.VERBOSE) return;
    if (size() != nodes->size()) {
        expect(false, "Mismatch between interpolator and mesh sizes: " + to_string(size()) + ", "
                + to_string(nodes->size()));
        return;
    }

    double rms_error = 0;
    int n_points = 0;
    for (int i = 0; i < nodes->size(); ++i)
        // look for the tetrahedral nodes whose potential == 0, i.e that are on the surface
        if ( nodes->get_marker(i) == TYPES.TETNODE && solutions[i].scalar == 0.0 &&
                c.inside_interesting_region((*nodes)[i]) )
        {
            double analyt = get_analyt_field(i).norm();
            double numerical = solutions[i].norm;

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

    const int n_atoms = solutions.size();
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

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
// surrounding hexahedral nodes
bool TetrahedronInterpolator::average_sharp_nodes() {
    vector<vector<unsigned int>> vorocells;
    mesh->calc_pseudo_3D_vorocells(vorocells);
    return LinearInterpolator<4>::average_sharp_nodes(vorocells, elems->stat.edgemax);
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
            append_solution(Solution(error_field));
    }

    // remove the spikes in the solution
    if (average_sharp_nodes())
        return true;

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : *elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
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
            append_solution(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : *elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
                return true;

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
            append_solution(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : *elems)
        for (int node : elem)
            if (solutions[node].scalar == error_field)
                return true;

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

// Precompute the data to tetrahedra to make later bcc calculations faster
void TetrahedronInterpolator::precompute() {
    const int n_elems = elems->size();
    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    reserve_precompute(n_elems);
    double d0, d1, d2, d3, d4;

    // Calculate tetrahedra neighbours
    for (int i = 0; i < n_elems; ++i)
        neighbours.push_back(elems->get_neighbours(i));

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroids.push_back(elems->get_centroid(i));

    // Store tetrahedra to become free to interpolate independently from mesh state
    for (int i = 0; i < n_elems; ++i)
        cells.push_back((*elems)[i]);

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
bool TetrahedronInterpolator::point_in_cell(const Vec3 &point, const int i) {
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

    // Ignore co-planar tetrahedra
    // no need to check because Tetgen guarantees non-co-planar tetrahedra
//    if (tet_not_valid[i]) return false;

    const Vec4 pt(point, 1);

    // If one of the barycentric coordinates is < zero, the point is outside the tetrahedron
    // Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
    if (det0[i] * pt.dotProduct(det1[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det2[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det3[i]) < zero) return false;
    if (det0[i] * pt.dotProduct(det4[i]) < zero) return false;

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

    return array<double,4> {bcc1, bcc2, bcc3, bcc4};
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
            if (nodes->get_marker(node) != 2)
                node_not_in_quads[node] = false;

    for (int node = 0; node < n_nodes; ++node)
        if (node_not_in_quads[node])
            solutions[node] = Solution(error_field);

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleFace face : *faces)
        for (int node : face)
            if (solutions[node].scalar == error_field)
                return true;

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
            append_solution(Solution(error_field));
    }

    // remove the spikes in the solution
    if (average_sharp_nodes())
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
bool TriangleInterpolator::average_sharp_nodes() {
    vector<vector<unsigned int>> vorocells;
    mesh->calc_pseudo_3D_vorocells(vorocells);
    return LinearInterpolator<3>::average_sharp_nodes(vorocells, mesh->faces.stat.edgemax);
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

    // Reserve memory for precomputation data
    reserve_precompute(n_faces);

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
            if (sface.neighbor((*faces)[j])) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
    }
}

// Calculate barycentric coordinates for a projection of a point inside the triangle
array<double,3> TriangleInterpolator::get_bcc(const Vec3& point, const int face) const {
    require(face >= 0 && face < vert0.size(), "Index out of bounds: " + to_string(face));

    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double u = tvec.dotProduct(pvec[face]);
    double v = qvec.dotProduct(norms[face]);

    return array<double,3> {1.0 - u - v, u, v};
}

// Check whether the projection of a point is inside the triangle
// It is separate routine from get_bcc to achieve better performance
bool TriangleInterpolator::point_in_cell(const Vec3& point, const int face) {
    require(face >= 0 && face < vert0.size(), "Index out of bounds: " + to_string(face));

    Vec3 tvec = point - vert0[face];
    double u = tvec.dotProduct(pvec[face]);
    if (u < zero || u > one) return false;     // Check first barycentric coordinate

    Vec3 qvec = tvec.crossProduct(edge1[face]);
    double v = qvec.dotProduct(norms[face]);
    if (v < zero || u + v > one) return false; // Check second & third barycentric coordinate

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
        if (scalar < error_field && nodes->get_marker(i) == TYPES.TETNODE)
            q_sum += scalar;
    }

    stringstream stream; stream << fixed << setprecision(3);
    stream << "Q / sum(" << scalar_label << ") = " << Q << " / " << q_sum << " = " << Q/q_sum;
    write_verbose_msg(stream.str());
}

template class LinearInterpolator<3> ;
template class LinearInterpolator<4> ;

} // namespace femocs
