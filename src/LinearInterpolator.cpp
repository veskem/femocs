/*
 * Interpolator.cpp
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#include "LinearInterpolator.h"
#include "Macros.h"

#include <float.h>

using namespace std;
namespace femocs {

// Initialize data vectors
LinearInterpolator::LinearInterpolator() {
    reserve(0);
    reserve_precompute(0);
}
;

// Print statistics about solution on node points
void LinearInterpolator::print_statistics() {
    if (!MODES.VERBOSE) return;

    const int n_atoms = size();
    Vec3 vec(0), rms_vec(0);
    double scalar = 0, rms_scalar = 0;
    int n_points = 0;

    for (int i = 0; i < n_atoms; ++i) {
        double s = solution[i].scalar;
        Vec3 v = solution[i].vector;
        if (s >= error_field) continue;

        vec += v; rms_vec += v * v;
        scalar += s; rms_scalar += s * s;
        n_points++;
    }

    vec *= (1.0 / n_points);
    rms_vec = Vec3(sqrt(rms_vec.x), sqrt(rms_vec.y), sqrt(rms_vec.z)) * (1.0 / n_points);
    scalar = scalar / n_points;
    rms_scalar = sqrt(rms_scalar) / n_points;

    cout << "  mean vector: \t" << vec << endl;
    cout << "   rms vector: \t" << rms_vec << endl;
    cout << "  mean & rms scalar: " << scalar << "\t" << rms_scalar << endl;
}

/* Return the mapping between tetrahedral & hexahedral mesh nodes,
 nodes & hexahedral elements and nodes & element's vertices  */
void LinearInterpolator::get_maps(dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh,
        vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs) {

    require(tria->n_vertices() > 0, "Invalid triangulation size!");
    const double eps = 1e-10;

    const int n_femocs_nodes = size();
    const int n_dealii_nodes = tria->n_used_vertices();
    const int n_verts_per_elem = dealii::GeometryInfo<3>::vertices_per_cell;

    vector<int> node2hex, node2vert;
    tet2hex.resize(n_femocs_nodes, -1);
    node2hex.resize(n_dealii_nodes);
    node2vert.resize(n_dealii_nodes);

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

// Force the solution on tetrahedral nodes to be the weighed average of the solutions on its voronoi cell nodes
void LinearInterpolator::average_tetnodes(const TetgenMesh &mesh) {
    vector<vector<unsigned int>> voro_cells = mesh.get_voronoi_cells();

    // loop through the tetrahedral nodes
    for (int i = 0; i < mesh.nodes.stat.n_tetnode; ++i) {
        if (solution[i].norm >= error_field) continue;

        Point3 tetnode = mesh.nodes[i];
        Vec3 vec(0);
        double w_sum = 0;

        // tetnode new solution will be the weighed average of the solutions on its voronoi cell nodes
        for (unsigned int v : voro_cells[i]) {
            double w = exp ( -1.0 * tetnode.distance(mesh.nodes[v]) );
            w_sum += w;
            vec += solution[v].vector * w;
        }

        solution[i].vector = vec * (1.0 / w_sum);
        solution[i].norm = solution[i].vector.norm();
    }
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool LinearInterpolator::extract_solution(fch::Laplace<3>* fem, const TetgenMesh &mesh) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = mesh.nodes.size();

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int i = 0; i < n_nodes; ++i)
        append( Atom(i, mesh.nodes[i], mesh.nodes.get_marker(i)) );

    // Precompute tetrahedra to make interpolation faster
    precompute_tetrahedra(mesh);

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tetNode2hexNode, cell_indxs, vert_indxs;
    get_maps(fem->get_triangulation(), fem->get_dof_handler(), tetNode2hexNode, cell_indxs, vert_indxs);

    vector<dealii::Tensor<1, 3>> ef = fem->get_efield(cell_indxs, vert_indxs); // get list of electric fields
    vector<double> pot = fem->get_potential(cell_indxs, vert_indxs); // get list of electric potentials

    require(ef.size() == pot.size(),
            "Mismatch of vector sizes: " + to_string(ef.size()) + ", " + to_string(pot.size()));

    int i = 0;
    for (int n : tetNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < ef.size(), "Invalid index: " + to_string(i));
            Vec3 elfield(ef[i][0], ef[i][1], ef[i][2]);
            elfield *= 0.1;   // convert V/nm  to  V/Angstom
            solution.push_back(Solution(elfield, 0.1*pot[i++]));
        }

        // In case of non-common node, store solution with error value
        else
            solution.push_back(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : mesh.elems)
        for (int node : elem)
            if (solution[node].scalar == error_field)
                return true;

    // force solution on tetrahedral nodes to be the weighed average of the solutions on its voronoi cell nodes
    average_tetnodes(mesh);

    return false;
}

// Extract the electric potential and electric field values on tetrahedral mesh nodes from FEM solution
bool LinearInterpolator::extract_solution(fch::CurrentsAndHeating<3>* fem, const TetgenMesh &mesh) {
    require(fem, "NULL pointer can't be handled!");
    const int n_nodes = mesh.nodes.size();

    // Copy the mesh nodes
    reserve(n_nodes);
    for (int n = 0; n < n_nodes; ++n)
        append( Atom(n, mesh.nodes[n], mesh.nodes.get_marker(n)) );

    // Precompute tetrahedra to make interpolation faster
    precompute_tetrahedra(mesh);

    // To make solution extraction faster, generate mapping between desired and available data sequences
    vector<int> tetNode2hexNode, cell_indxs, vert_indxs;
    get_maps(fem->get_triangulation(), fem->get_dof_handler(), tetNode2hexNode, cell_indxs, vert_indxs);

    vector<dealii::Tensor<1, 3>> rho = fem->get_current(cell_indxs, vert_indxs); // get list of current densities
    vector<double> temperature = fem->get_temperature(cell_indxs, vert_indxs); // get list of temperatures

    require(rho.size() == temperature.size(),
            "Mismatch of vector sizes: " + to_string(rho.size()) + ", "
                    + to_string(temperature.size()));

    int i = 0;
    for (int n : tetNode2hexNode) {
        // If there is a common node between tet and hex meshes, store actual solution
        if (n >= 0) {
            require(i < rho.size(), "Invalid index: " + to_string(i));
            Vec3 current(rho[i][0], rho[i][1], rho[i][2]);
            solution.push_back(Solution(current, temperature[i++]));
        }

        // In case of non-common node, store solution with error value
        else
            solution.push_back(Solution(error_field));
    }

    // Check for the error values in the mesh nodes
    // Normally there should be no nodes in the mesh elements that have the error value
    for (SimpleElement elem : mesh.elems)
        for (int node : elem)
            if (solution[node].scalar == error_field)
                return true;

    return false;
}

// Reserve memory for interpolation data
void LinearInterpolator::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + to_string(N));
    atoms.clear();
    solution.clear();

    atoms.reserve(N);
    solution.reserve(N);
}

// Reserve memory for pre-compute data
void LinearInterpolator::reserve_precompute(const int N) {
    tetrahedra.clear();
    tetneighbours.clear();
    centroid.clear();
    det0.clear();
    det1.clear();
    det2.clear();
    det3.clear();
    det4.clear();
    tet_not_valid.clear();

    tetrahedra.reserve(N);
    tetneighbours.reserve(N);
    centroid.reserve(N);
    det0.reserve(N);
    det1.reserve(N);
    det2.reserve(N);
    det3.reserve(N);
    det4.reserve(N);
    tet_not_valid.reserve(N);
}

// Precompute the data to tetrahedra to make later bcc calculations faster
void LinearInterpolator::precompute_tetrahedra(const TetgenMesh &mesh) {
    const int n_elems = mesh.elems.size();
    double d0, d1, d2, d3, d4;

    reserve_precompute(n_elems);

    // Copy tetrahedra
    for (int i = 0; i < n_elems; ++i)
        tetrahedra.push_back(mesh.elems[i]);

    // Copy tetrahedra neighbours
    for (int i = 0; i < n_elems; ++i)
        tetneighbours.push_back(mesh.elems.get_neighbours(i));

    // Calculate centroids of elements
    for (int i = 0; i < n_elems; ++i)
        centroid.push_back(mesh.elems.get_centroid(i));

    /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
     * barycentric coordinate of tetrahedra using the relations below */
    for (SimpleElement se : tetrahedra) {
        Vec3 v1 = mesh.nodes.get_vec(se[0]);
        Vec3 v2 = mesh.nodes.get_vec(se[1]);
        Vec3 v3 = mesh.nodes.get_vec(se[2]);
        Vec3 v4 = mesh.nodes.get_vec(se[3]);

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
bool LinearInterpolator::point_in_tetrahedron(const Point3 &point, const int i) {
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

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
Vec4 LinearInterpolator::get_bcc(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < det0.size(), "Index out of bounds: " + to_string(elem));

    Vec4 pt(point, 1);
    double bcc1 = det0[elem] * pt.dotProduct(det1[elem]);
    double bcc2 = det0[elem] * pt.dotProduct(det2[elem]);
    double bcc3 = det0[elem] * pt.dotProduct(det3[elem]);
    double bcc4 = det0[elem] * pt.dotProduct(det4[elem]);

    return Vec4(bcc1, bcc2, bcc3, bcc4);
}

// Find the element which contains the point or is the closest to it
int LinearInterpolator::locate_element(const Point3 &point, const int elem_guess) {
    const int n_elems = det0.size();

    // Check the guessed element
    if (point_in_tetrahedron(point, elem_guess)) return elem_guess;

    // Check the 1st nearest neighbours of guessed element
    for (int nbor : tetneighbours[elem_guess])
        if (nbor >= 0 && point_in_tetrahedron(point, nbor)) return nbor;

    // Check the 2nd, 3rd & 4th nearest neighbours of guessed element
    // going further than 4th neighbour starts slowing things down

    // Mark the neighbour ranks of elements wrt to guessed element
    vector<int> nbor_rank(n_elems);
    for (int nbor1 : tetneighbours[elem_guess]) {
        if (nbor1 < 0) continue;
        for (int nbor2 : tetneighbours[nbor1]) {
            if (nbor2 < 0) continue;
            for (int nbor3 : tetneighbours[nbor2]) {
                if (nbor3 < 0) continue;
                for (int nbor4 : tetneighbours[nbor3]) {
                    if (nbor4 < 0) continue;
                    nbor_rank[nbor4] = 4;
                }
                nbor_rank[nbor3] = 3;
            }
            nbor_rank[nbor2] = 2;
        }
        nbor_rank[nbor1] = 1;
    }

    // Perform the check on the neighbours
    // checking ranks separately doesn't give any benefit in speed
    for (int elem = 0; elem < n_elems; ++elem)
        if (nbor_rank[elem] > 1 && point_in_tetrahedron(point, elem)) return elem;

    // If no success, loop through all the elements
    double min_distance2 = DBL_MAX;
    int min_index = 0;

    for (int elem = 0; elem < n_elems; ++elem) {
        // If correct element is found, we're done
        if (point_in_tetrahedron(point, elem))
            return elem;

        // Otherwise look for the element whose centroid is closest to the point
        else {
            double distance2 = point.distance2(centroid[elem]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = elem;
            }
        }
    }

    // If no perfect element found, return the best
    // indicate the imperfectness with the minus sign
    return -min_index;
}

// Calculate interpolation or all the data for point inside or near the elem-th tetrahedron
Solution LinearInterpolator::get_solution(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < tetrahedra.size(), "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);

    SimpleElement selem = tetrahedra[elem];

    // Interpolate electric field
    Vec3 elfield_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        elfield_i += solution[selem[i]].vector * bcc[i];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        potential_i += solution[selem[i]].scalar * bcc[i];

    return Solution(elfield_i, potential_i);
}

// Return full solution on i-th node
Solution LinearInterpolator::get_solution(const int i) {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solution[i];
}

// Calculate interpolation for vector data for point inside or near the elem-th tetrahedron
Vec3 LinearInterpolator::get_vector(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < tetrahedra.size(), "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);
    SimpleElement selem = tetrahedra[elem];

    // Interpolate electric field
    Vec3 vec(0.0);
    for (int i = 0; i < selem.size(); ++i)
        vec += solution[selem[i]].vector * bcc[i];

    return vec;
}

// Return vector component of solution on i-th node
Vec3 LinearInterpolator::get_vector(const int i) {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solution[i].vector;
}

// Calculate interpolation for scalar data for point inside or near the elem-th tetrahedron
double LinearInterpolator::get_scalar(const Point3 &point, const int elem) {
    require(elem >= 0 && elem < tetrahedra.size(), "Index out of bounds: " + to_string(elem));

    // Get barycentric coordinates of point in tetrahedron
    Vec4 bcc = get_bcc(point, elem);
    SimpleElement selem = tetrahedra[elem];

    // Interpolate potential
    double potential_i(0.0);
    for (int i = 0; i < selem.size(); ++i)
        potential_i += solution[selem[i]].scalar * bcc[i];

    return potential_i;
}

// Return scalar component of solution on i-th node
double LinearInterpolator::get_scalar(const int i) {
    require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
    return solution[i].scalar;
}

// Compile data string from the data vectors for file output
string LinearInterpolator::get_data_string(const int i) const {
    if (i < 0)
        return "LinearInterpolator properties=id:R:1:pos:R:3:marker:R:1:force:R:3:enorm:R:1:potential:R:1";

    ostringstream strs;
    strs << atoms[i] << " " << solution[i];
    return strs.str();
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

} // namespace femocs
