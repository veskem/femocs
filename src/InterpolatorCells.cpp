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
        mesh(NULL), norm_label("vector_norm"), scalar_label("scalar") { reserve(0); }

InterpolatorNodes::InterpolatorNodes(const string &nl, const string &sl) :
        mesh(NULL), norm_label(nl), scalar_label(sl) { reserve(0); }

void InterpolatorNodes::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + to_string(N));

    markers.clear();
    markers.reserve(N);
    vertices.clear();
    vertices.reserve(N);
    solutions.clear();
    solutions.reserve(N);
}

void InterpolatorNodes::precompute() {
    const int n_nodes = mesh->nodes.size();
    reserve(n_nodes);

    for (int i = 0; i < n_nodes; ++i) {
        // Store the coordinates of mesh nodes
        vertices.push_back(mesh->nodes[i]);
        // Store node markers
        markers.push_back(mesh->nodes.get_marker(i));
    }
}

void InterpolatorNodes::write(const string &file_name) const {
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

void InterpolatorNodes::write_xyz(ofstream& out) const {
    const int n_nodes = vertices.size();
    expect(n_nodes, "Zero nodes detected!");

    out << n_nodes << endl;
    out << "Interpolator properties=id:I:1:pos:R:3:marker:I:1:" <<
            "force:R:3:" << norm_label << ":R:1:" << scalar_label << ":R:1" << endl;

    for (int i = 0; i < n_nodes; ++i)
        out << i << " " << vertices[i] << " " << markers[i] << " " << solutions[i] << endl;
}

void InterpolatorNodes::write_vtk(ofstream& out) const {
    const int n_nodes = size();
    const int n_cells = size();
    const int dim = 1;
    const int celltype = get_cell_type();

    expect(n_nodes, "Zero nodes detected!");

    out << "# vtk DataFile Version 3.0\n";
    out << "# InterpolatorNodes data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << vertices[i] << "\n";

    // Output # vertices and vertex indices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << i << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";

    write_point_data(out);
}

void InterpolatorNodes::write_point_data(ofstream& out) const {
    const int n_nodes = size();
    out << "\nPOINT_DATA " << n_nodes << "\n";

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

    out << "VECTORS " << norm_label << "_vec double\n";
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
    require(N >= 0, "Invalid number of points: " + to_string(N));

    centroids.clear();
    centroids.reserve(N);
    cells.clear();
    cells.reserve(N);
    markers = vector<int>(N);
    neighbours = vector<vector<int>>(N);
}

template<int dim>
void InterpolatorCells<dim>::write(const string &file_name) const {
    if (!MODES.WRITEFILE) return;

    ofstream outfile;
    outfile.setf(std::ios::scientific);
    outfile.precision(6);

    string ftype = get_file_type(file_name);
    outfile.open(file_name);
    require(outfile.is_open(), "Can't open a file " + file_name);

    if (ftype == "vtk")
        write_vtk(outfile);
    else
        require(false, "Unsupported file type: " + ftype);

    outfile.close();
}

template<int dim>
void InterpolatorCells<dim>::write_vtk(ofstream& out) const {
    const int n_nodes = nodes->size(); // nodes->stat.n_tetnode;
    expect(n_nodes, "Zero nodes detected!");

    const int celltype = get_cell_type();
    const int n_cells = cells.size();

    out << "# vtk DataFile Version 3.0\n";
    out << "# InterpolatorCells data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // Output the point coordinates
    out << "POINTS " << n_nodes << " double\n";
    for (int i = 0; i < n_nodes; ++i)
        out << nodes->get_vertex(i) << "\n";

    // Output # vertices and vertex indices
    out << "\nCELLS " << n_cells << " " << (1+dim) * n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << dim << " " << cells[i] << "\n";

    // Output cell types
    out << "\nCELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i)
        out << celltype << "\n";

    // Associate solution with the nodes
    nodes->write_point_data(out);

    // Write data associated with the cells
    write_cell_data(out);
}

template<int dim>
void InterpolatorCells<dim>::write_cell_data(ofstream& out) const {
    const int n_cells = cells.size();
    out << "\nCELL_DATA " << n_cells << "\n";

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
    // Check the guessed cell
    Vec3 vec_point(point);
    if (point_in_cell(vec_point, cell_guess)) return cell_guess;

    vector<bool> cell_checked = vector_not(&markers, 0);
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

template<int dim>
void InterpolatorCells<dim>::get_weights(array<double,dim>& weights, const Point3 &point, const SimpleCell<dim>& scell) const {
    double w_sum = 0;
    for (int i = 0; i < dim; ++i) {
        double w = exp( decay_factor * point.distance(nodes->get_vertex(scell[i])) );
        weights[i] = w;
        w_sum += w;
    }

    // check for zero and nan
    require(w_sum > 0 && w_sum == w_sum, "Invalid interpolation weight: " + to_string(w_sum));
    w_sum = 1.0 / w_sum;
    for (int i = 0; i < dim; ++i)
        weights[i] *= w_sum;
}

template<int dim>
Solution InterpolatorCells<dim>::interp_solution(const Point3 &point, const int c) const {
    const int cell = abs(c);
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

    SimpleCell<dim> scell = cells[cell];

    // calculate weights or barycentric coordinates
    array<double,dim> weights;
    if (c >= 0) get_shape_functions(weights, Vec3(point), cell);
    else get_weights(weights, point, scell);

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
    require(cell < cells.size(), "Index out of bounds: " + to_string(cell));

    // calculate shape functions and their gradients
    array<double, dim> sf;
    array<Vec3, dim> sfg;
    get_shape_functions(sf, point, cell);
    get_shape_fun_grads(sfg, point, cell);

    // using them as weights, interpolate scalar data and its gradient
    SimpleCell<dim> scell = cells[cell];
    Vec3 vector_i(0.0);
    double scalar_i(0.0);

    for (int i = 0; i < dim; ++i) {
        double scalar = nodes->get_scalar(scell[i]);
        vector_i += sfg[i] * scalar;
        scalar_i += sf[i] * scalar;
    }

    return Solution(vector_i, scalar_i);
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
 *  ======================= LinearTriangles ========================
 * ================================================================== */

LinearTriangles::LinearTriangles() :
        InterpolatorCells<3>(), tris(NULL) {}

LinearTriangles::LinearTriangles(const InterpolatorNodes* n) :
        InterpolatorCells<3>(n), tris(NULL) {}

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
    require(mesh && tris, "NULL pointers can't be used!");
    const int n_faces = tris->size();
    const int n_nodes = nodes->size();

    expect(n_faces > 0, "Interpolator expects non-empty mesh!");

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Store the constant for smoothing
    decay_factor = -1.0 / tris->stat.edgemax;

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        SimpleFace sface = (*tris)[i];

        Vec3 v0 = mesh->nodes.get_vec(sface[0]);
        Vec3 v1 = mesh->nodes.get_vec(sface[1]);
        Vec3 v2 = mesh->nodes.get_vec(sface[2]);

        Vec3 e1 = v1 - v0;   // edge1 of triangle
        Vec3 e2 = v2 - v0;   // edge2 of triangle
        Vec3 pv = tris->get_norm(i).crossProduct(e2);
        double i_det = 1.0 / e1.dotProduct(pv);

        vert0.push_back(v0);
        edge1.push_back(e1 * i_det);
        edge2.push_back(e2);
        pvec.push_back(pv * i_det);

        // store triangles to get rid of dependence on mesh state
        cells.push_back(sface);
        // calculate norms of triangles
        norms.push_back(tris->get_norm(i));
        // store max distance from given triangle
        max_distance.push_back(e2.norm());
        // calculate centroids of triangles
        centroids.push_back(tris->get_centroid(i));

        // calculate the neighbour list for triangles
        for (int j = i+1; j < n_faces; ++j)
            if (sface.edge_neighbor((*tris)[j])) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
    }
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
        get_shape_functions(weights, Vec3(point), cell);
        // store the cell to make next step faster
        atom2cell[i] = cell;
        // append barycentric weight to the sum of all weights from given node
        int j = 0;
        for (int node : cells[cell])
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
        for (int node : cells[cell])
            scalars[i] += nodes->get_scalar(node) * weights[j++] / bcc_sum[node];
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

void LinearTriangles::get_shape_functions(array<double,3>& sf, const Vec3& point, const int face) const {
    Vec3 tvec = point - vert0[face];
    Vec3 qvec = tvec.crossProduct(edge1[face]);
    const double v = tvec.dotProduct(pvec[face]);
    const double w = qvec.dotProduct(norms[face]);
    const double u = 1.0 - v - w;
    sf = {zero + u, zero + v, zero + w};
}

int LinearTriangles::near_surface(const Vec3& point, const double r_cut) const {
    require(r_cut > 0, "Invalid distance from surface: " + to_string(r_cut));

    for (int face = 0; face < cells.size(); ++face) {
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

void LinearTriangles::write_cell_data(ofstream& out) const {
    InterpolatorCells<3>::write_cell_data(out);

    // write face norms
    out << "VECTORS norm double\n";
    for (int i = 0; i < cells.size(); ++i)
        out << norms[i] << "\n";
}

/* ==================================================================
 *  ====================== QuadraticTriangles ======================
 * ================================================================== */

QuadraticTriangles::QuadraticTriangles() :
        InterpolatorCells<6>(), tris(NULL), lintri(NULL) {}

QuadraticTriangles::QuadraticTriangles(const InterpolatorNodes* n, const LinearTriangles* l) :
    InterpolatorCells<6>(n), tris(NULL), lintri(l) {}

void QuadraticTriangles::reserve(const int N) {
    InterpolatorCells<6>::reserve(N);
}

void QuadraticTriangles::precompute() {
    require(mesh && tris && lintri, "NULL pointers can't be used!");
    const int n_faces = tris->size();
    require(n_faces == lintri->size(), "QuadraticTriangles requires LinearTriangles to be pre-computed!");

    // Reserve memory for precomputation data
    reserve(n_faces);

    // Store the constant for smoothing
    decay_factor = -1.0 / tris->stat.edgemax;

    // Loop through all the faces
    for (int i = 0; i < n_faces; ++i) {
        // store triangles to get rid of dependence on mesh state
        cells.push_back(get_cell(i));

        // calculate and store centroids of triangles
        centroids.push_back(tris->get_centroid(i));

        // calculate the neighbour list for triangles
        SimpleFace sface = (*tris)[i];
        for (int j = i + 1; j < n_faces; ++j)
            if (sface.edge_neighbor((*tris)[j])) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
    }
}

bool QuadraticTriangles::point_in_cell(const Vec3& point, const int face) const {
    return lintri->point_in_cell(point, face);
}

void QuadraticTriangles::get_shape_functions(array<double,6>& sf, const Vec3& point, const int face) const {
    array<double,3> bcc;
    lintri->get_shape_functions(bcc, point, face);

    const double u = bcc[0];
    const double v = bcc[1];
    const double w = bcc[2];

    sf[0] = u * (2 * u - 1);
    sf[1] = v * (2 * v - 1);
    sf[2] = w * (2 * w - 1);
    sf[3] = 4 * u * v;
    sf[4] = 4 * v * w;
    sf[5] = 4 * w * u;

//    sf = {u, v, w, 0, 0, 0};
}

SimpleCell<6> QuadraticTriangles::get_cell(const int tri) const {
    require(tri >= 0 && tri < tris->size(), "Invalid index: " + to_string(tri));
    if (mesh->quads.size() == 0)
        return QuadraticTri(0);

    const int n_quads_per_tri = 3;
    array<vector<unsigned>,3> edge_nodes;

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
    const int n_nodes = this->nodes->size();

    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    double d0, d1, d2, d3, d4;

    reserve(n_elems);

    // Store the constant for smoothing
    decay_factor = -1.0 / tets->stat.edgemax;

    for (int i = 0; i < n_elems; ++i) {
        SimpleElement se = (*tets)[i];

        // Calculate tetrahedra neighbours
        neighbours[i] = tets->get_neighbours(i);
        // Calculate centroids of tetrahedra
        centroids.push_back(tets->get_centroid(i));
        // Store tetrahedra to get rid of mesh dependency
        cells.push_back(se);

        /* Calculate main and minor determinants for 1st, 2nd, 3rd and 4th
         * barycentric coordinate of tetrahedra using the relations below */

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
    require(i >= 0 && i < det0.size(), "Index out of bounds: " + to_string(i));

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

void LinearTetrahedra::get_shape_functions(array<double,4>& sf, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < det0.size(), "Index out of bounds: " + to_string(tet));

    const Vec4 pt(point, 1);
    sf[0] = zero + det0[tet] * pt.dotProduct(det1[tet]);
    sf[1] = zero + det0[tet] * pt.dotProduct(det2[tet]);
    sf[2] = zero + det0[tet] * pt.dotProduct(det3[tet]);
    sf[3] = zero + det0[tet] * pt.dotProduct(det4[tet]);
}

void LinearTetrahedra::get_shape_fun_grads(array<Vec3, 4>& sfg, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < cells.size(), "Index out of bounds: " + to_string(tet));

    // The gradient of shape function inside linear tetrahedron does not depend
    // on the point of interest, therefore its value could be precalculated and
    // used upon request. However, for the consistency with other cells it's not
    // done now. Maybe in the future, if someone wants to use this function very heavily.

    // calculate 4x4 Jacobian
    // the upmost row consists of 1-s and is not used, therefore no need to store it
    SimpleCell<4> cell = cells[tet];
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
    sfg = {Jinv[0], Jinv[1], Jinv[2], Jinv[3]};
}

void LinearTetrahedra::narrow_search_to(const int region) {
    const int n_cells = cells.size();
    require(n_cells == tets->size(), "LinearTetrahedra must be intialized before narrowing search!");

    const int surf_end = mesh->nodes.indxs.surf_end;

    if (region == TYPES.VACUUM) {
        for (int i = 0; i < n_cells; ++i)
            markers[i] = tets->get_marker(i) != TYPES.VACUUM;
    } else if (region == TYPES.BULK) {
        for (int i = 0; i < n_cells; ++i)
            markers[i] = tets->get_marker(i) == TYPES.VACUUM;
    } else if (region == TYPES.SURFACE) {
        for (int i = 0; i < n_cells; ++i) {
            SimpleElement se = cells[i];
            markers[i] = se[0] > surf_end && se[1] > surf_end && se[2] > surf_end && se[3] > surf_end;
        }
    } else
        require(false, "Unimplemented region: " + to_string(region));
}

/* ==================================================================
 *  ===================== QuadraticTetrahedra ======================
 * ================================================================== */

QuadraticTetrahedra::QuadraticTetrahedra() :
        InterpolatorCells<10>(), tets(NULL), lintet(NULL) {}

QuadraticTetrahedra::QuadraticTetrahedra(const InterpolatorNodes* n, const LinearTetrahedra* l) :
    InterpolatorCells<10>(n), tets(NULL), lintet(l) {}

void QuadraticTetrahedra::reserve(const int N) {
    InterpolatorCells<10>::reserve(N);
}

void QuadraticTetrahedra::precompute() {
    require(mesh && tets && lintet, "NULL pointers can't be used!");

    const int n_elems = tets->size();
    require(n_elems == lintet->size(), "QuadraticTetrahedra requires LinearTetrahedra to be pre-computed!");

    reserve(n_elems);

    // Store the constant for smoothing
    this->decay_factor = -1.0 / tets->stat.edgemax;

    // Loop through all the tetrahedra
    for (int i = 0; i < n_elems; ++i) {
        // Calculate tetrahedra neighbours
        neighbours[i] = tets->get_neighbours(i);
        // Calculate centroids of tetrahedra
        centroids.push_back(tets->get_centroid(i));
        // Calculate and store 10-noded tetrahedra
        cells.push_back(get_cell(i));
    }
}

bool QuadraticTetrahedra::point_in_cell(const Vec3& point, const int face) const {
    return lintet->point_in_cell(point, face);
}

void QuadraticTetrahedra::get_shape_functions(array<double,10>& sf, const Vec3& point, const int tet) const {
    array<double,4> bcc;
    lintet->get_shape_functions(bcc, point, tet);

    const double b1 = bcc[0];
    const double b2 = bcc[1];
    const double b3 = bcc[2];
    const double b4 = bcc[3];

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
void QuadraticTetrahedra::get_shape_fun_grads(array<Vec3, 10>& sfg, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < cells.size(), "Index out of bounds: " + to_string(tet));

    array<double,4> bcc;
    lintet->get_shape_functions(bcc, point, tet);

    // Store tetrahedron nodes for more convenient access
    SimpleCell<10> cell = cells[tet];
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
    sfg = {
        Jinv[0] * (4*bcc[0]-1),
        Jinv[1] * (4*bcc[1]-1),
        Jinv[2] * (4*bcc[2]-1),
        Jinv[3] * (4*bcc[3]-1),
        Jinv[0]*bcc[1] + Jinv[1]*bcc[0],
        Jinv[1]*bcc[2] + Jinv[2]*bcc[1],
        Jinv[0]*bcc[2] + Jinv[2]*bcc[0],
        Jinv[0]*bcc[3] + Jinv[3]*bcc[0],
        Jinv[1]*bcc[3] + Jinv[3]*bcc[1],
        Jinv[2]*bcc[3] + Jinv[3]*bcc[2]
    };
    for (int i = 4; i < 10; ++i)
        sfg[i] *= 4.0;
}

void QuadraticTetrahedra::get_shape_fun_grads_slow(array<Vec3, 10>& sfg, const Vec3& point, const int tet) const {
    require(tet >= 0 && tet < cells.size(), "Index out of bounds: " + to_string(tet));

    array<double,4> bcc;
    lintet->get_shape_functions(bcc, point, tet);

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

    // calculate gradient of shape functions in xyz-space
    for (int k = 0; k < 10; ++k) {
        sfg[k] = Vec3(0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j)
                sfg[k][i] += Jinv[j][i+1] * dN[k][j];
    }
}

void QuadraticTetrahedra::test_shape_funs() {
    array<Vec3, 10> sfg1, sfg2;

    cout << fixed << setprecision(3);

    int tet = 0;
    QuadraticTet cell = cells[tet];
    for (int i = 0; i < 4; ++i) {
        cout << "\ni = " << i << endl;
        Point3 point = mesh->nodes[cell[i]];
        get_shape_fun_grads(sfg1, point, tet);
        get_shape_fun_grads_slow(sfg2, point, tet);

        cout << "delta sfg = " << endl;
        for (int i = 0; i < 10; ++i)
            cout << sfg1[i] - sfg2[i] << endl;
    }
}

SimpleCell<10> QuadraticTetrahedra::get_cell(const int tet) const {
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
    require(N >= 0, "Invalid number of points: " + to_string(N));

    markers = vector<int>(N);
    cells.clear(); cells.reserve(N);
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
    require(mesh->tets.size() == lintet->size(), "LinearHexahedra requires LinearTetrahedra to be pre-computed!");

    const int n_elems = hexs->size();
    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    reserve(n_elems);

    decay_factor = -1.0 / mesh->tets.stat.edgemax;

    // Loop through all the hexahedra
    for (int i = 0; i < n_elems; ++i) {
        SimpleHex cell = (*hexs)[i];

        // store hexahedra
        cells.push_back(cell);

        // pre-calculate data to make iterpolation faster
        const Vec3 x1 = mesh->nodes.get_vec(cell[0]);
        const Vec3 x2 = mesh->nodes.get_vec(cell[1]);
        const Vec3 x3 = mesh->nodes.get_vec(cell[2]);
        const Vec3 x4 = mesh->nodes.get_vec(cell[3]);
        const Vec3 x5 = mesh->nodes.get_vec(cell[4]);
        const Vec3 x6 = mesh->nodes.get_vec(cell[5]);
        const Vec3 x7 = mesh->nodes.get_vec(cell[6]);
        const Vec3 x8 = mesh->nodes.get_vec(cell[7]);

        f0s.push_back( Vec3((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8) / 8.0) );
        f1s.push_back( Vec3(((x1*-1) + x2 + x3 - x4 - x5 + x6 + x7 - x8) / 8.0) );
        f2s.push_back( Vec3(((x1*-1) - x2 + x3 + x4 - x5 - x6 + x7 + x8) / 8.0) );
        f3s.push_back( Vec3(((x1*-1) - x2 - x3 - x4 + x5 + x6 + x7 + x8) / 8.0) );
        f4s.push_back( Vec3((x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8) / 8.0) );
        f5s.push_back( Vec3((x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8) / 8.0) );
        f6s.push_back( Vec3((x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8) / 8.0) );
        f7s.push_back( Vec3(((x1*-1) + x2 - x3 + x4 + x5 - x6 + x7 - x8) / 8.0) );
    }

    // make the markers to correspond to lintet
    for (int i = 0; i < n_elems; ++i)
        markers[i] = lintet->get_marker(hexs->to_tet(i));

    // store the mapping between femocs and deal.ii hexahedra
    int deal_hex_index = 0;
    map_femocs2deal = vector<int>(n_elems);
    for (int i = 0; i < n_elems; ++i) {
        if (hexs->get_marker(i) > 0)
            map_femocs2deal[i] = deal_hex_index++;
        else
            map_femocs2deal[i] = -1;
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
        require(D != 0, "Invalid determinant: " + to_string(D));

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

void LinearHexahedra::get_shape_functions(array<double,8>& sf, const Vec3& point, const int hex) const {
    require(hex >= 0 && hex < cells.size(), "Index out of bounds: " + to_string(hex));

    // Before calculating the shape function,
    // map the point from Cartesian xyz-coordinates to natural uvw-coordinates.
    double u, v, w;
    project_to_nat_coords(u, v, w, point, hex);

    // use natural coordinates to calculate shape functions
    sf = {
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

/* For more information, see the lecture materials in
 * https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
 */
void LinearHexahedra::get_shape_fun_grads(array<Vec3, 8>& sfg, const Vec3& point, const int hex) const {
    require(hex >= 0 && hex < cells.size(), "Index out of bounds: " + to_string(hex));
    double u, v, w;
    project_to_nat_coords(u, v, w, point, hex);

    // Transfer the coordinates of hex nodes into a matrix
    SimpleHex shex = mesh->hexs[hex];
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

    // calculate gradient of shape functions in xyz-space
    for (int k = 0; k < 8; ++k) {
        sfg[k] = Vec3(0);
        for (int i = 0; i < 3; ++i)
            sfg[k] += Jinv[i] * dN[k][i];
    }
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
void LinearHexahedra::get_shape_fun_grads(array<Vec3, 8>& sfg, const int hex, const int node) const {
    require(hex >= 0 && hex < size(), "Invalid hex: " + to_string(hex));
    require(node >= 0 && node < 8, "Invalid node: " + to_string(node));

    double u, v, w; // natural coordinates
    int n1, n2, n3; // indices of nearest neighbouring nodes of the node
    project_to_nat_coords(u, v, w, n1, n2, n3, node);

    Vec3 uvw(u, v, w);
    int nnn[3] = {n1, n2, n3};
    SimpleHex shex = mesh->hexs[hex];

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
    sfg = array<Vec3,8>();
    for (int i = 0; i < 3; ++i) {
        sfg[node][i] = 0.5 * (Jinv[i].dotProduct(uvw));
        for (int j = 0; j < 3; ++j)
            sfg[nnn[j]][i] = -0.5 * Jinv[i][j] * uvw[j];
    }
}

void LinearHexahedra::test_shape_funs() {
    array<Vec3, 8> sfg1, sfg2;

    int hex = 0;
    SimpleHex shex = mesh->hexs[hex];

    cout << fixed << setprecision(3);
    for (int i = 0; i < 1; ++i) {
        Point3 point = mesh->nodes[shex[i]];
        get_shape_fun_grads(sfg1, point, hex);

        cout << endl << i << endl;
        for (int j = 0; j < 8; ++j)
            cout << sfg2[j] - sfg1[j] << endl;
    }
}

/*
 * Function uses the fact that hexahedra are always uniquely tied with the nodes of tetrahedra.
 * Therefore, knowing the barycentric coordinates inside a tetrahedron
 * makes it possible to use them to determine, in which part of the tetrahedron and therefore
 * also inside which hexahedron the point is located.
 */
int LinearHexahedra::locate_cell(const Point3 &point, const int cell_guess) const {
    static constexpr int n_hexs_per_tet = 4;

    int tet_index = abs(cell_guess / n_hexs_per_tet);
    tet_index = lintet->locate_cell(point, tet_index);
    int sign = 1;
    if (tet_index < 0) sign = -1;
    tet_index = abs(tet_index);

    // calculate barycentric coordinates for a point
    array<double,4> bcc;
    lintet->get_shape_functions(bcc, point, tet_index);

    // all the ratios below are == 1, if the point is exactly on the boundary of two hexahedra
    // and 0 or inf, if point is on the face of a tetrahedron.
    // if ratio is > 1, the point is closer to 1st corresponding node, if < 1, it's closer to 2nd node
    const double b1b2 = bcc[0] / bcc[1];
    const double b1b3 = bcc[0] / bcc[2];
    const double b1b4 = bcc[0] / bcc[3];
    const double b2b3 = bcc[1] / bcc[2];
    const double b2b4 = bcc[1] / bcc[3];
    const double b3b4 = bcc[2] / bcc[3];

    // point inside a hex connected to 1st tetrahedral node ?
    if (b1b2 >= 1 && b1b3 >= 1 && b1b4 >= 1)
        return sign * (n_hexs_per_tet * tet_index + 0);
    // point inside a hex connected to 2nd tetrahedral node ?
    if (b1b2 <= 1 && b2b3 >= 1 && b2b4 >= 1)
        return sign * (n_hexs_per_tet * tet_index + 1);
    // point inside a hex connected to 3rd tetrahedral node ?
    if (b1b3 <= 1 && b2b3 <= 1 && b3b4 >= 1)
        return sign * (n_hexs_per_tet * tet_index + 2);
    // point inside a hex connected to 4th tetrahedral node ?
    if (b1b4 <= 1 && b2b4 <= 1 && b3b4 <= 1)
        return sign * (n_hexs_per_tet * tet_index + 3);

    return -1;
}

/* ==================================================================
 *  ====================== LinearQuadrangles =======================
 * ================================================================== */

LinearQuadrangles::LinearQuadrangles() :
        InterpolatorCells<4>(), quads(NULL), lintri(NULL) {}

LinearQuadrangles::LinearQuadrangles(const InterpolatorNodes* n, const LinearTriangles* l) :
        InterpolatorCells<4>(n), quads(NULL), lintri(l) {}

void LinearQuadrangles::reserve(const int N) {
    require(N >= 0, "Invalid number of points: " + to_string(N));

    markers = vector<int>(N);
    cells.clear(); cells.reserve(N);
}

void LinearQuadrangles::precompute() {
    require(mesh && quads && lintri, "NULL pointers can't be used!");
    require(mesh->tris.size() == lintri->size(), "LinearQuadrangles requires LinearTriangles to be pre-computed!");

    const int n_elems = quads->size();
    expect(n_elems > 0, "Interpolator expects non-empty mesh!");
    reserve(n_elems);

    // Loop through all the hexahedra
    for (int i = 0; i < n_elems; ++i) {
        SimpleQuad cell = (*quads)[i];

        // store hexahedra
        cells.push_back(cell);
    }

    // make the markers to correspond to lintri
    for (int i = 0; i < n_elems; ++i)
        markers[i] = lintri->get_marker(quads->to_tri(i));
}

Solution LinearQuadrangles::interp_solution(const Point3 &point, const int c) const {
    int cell;
    if (c >= 0) cell = quads->to_tri(c);
    else cell = -quads->to_tri(abs(c));

    return lintri->interp_solution(point, cell);
}

/*
 * The code for mapping the point from Cartesian to natural coordinate space was taken from
 * https://www.gamedev.net/forums/topic/596392-uv-coordinate-on-a-2d-quadrilateral/
 */
void LinearQuadrangles::get_shape_functions(array<double,4>& sf, const Vec3& p, const int quad) const {
    require(false, "LinearQuadrangles::get_shape_functions is not implemented!");

    require(quad >= 0 && quad < cells.size(), "Index out of bounds: " + to_string(quad));

    // Before calculating the shape function,
    // map the point from Cartesian xyz-coordinates to natural uv-coordinates.
    // In natural coordinate system, each coordinate is within limits [-1, 1].

    SimpleQuad squad = cells[quad];
    Vec3 a = mesh->nodes.get_vec(squad[0]);
    Vec3 b = mesh->nodes.get_vec(squad[1]);
    Vec3 c = mesh->nodes.get_vec(squad[2]);
    Vec3 d = mesh->nodes.get_vec(squad[3]);

    double C = (a.y - p.y) * (d.x - p.x) - (a.x - p.x) * (d.y - p.y);
    double B = (a.y - p.y) * (c.x - d.x) + (b.y - a.y) * (d.x - p.x) - (a.x - p.x) * (c.y - d.y) - (b.x - a.x) * (d.y - p.y);
    double A = (b.y - a.y) * (c.x - d.x) - (b.x - a.x) * (c.y - d.y);
    double D = B * B - 4 * A * C;

    double u = (-B - sqrt(D)) / (2 * A);
    double p1x = a.x + (b.x - a.x) * u;
    double p2x = d.x + (c.x - d.x) * u;
    double v = (p.x - p1x) / (p2x - p1x);

    // use natural coordinates to calculate shape functions
//    sf = {
//            (1 - u) * (1 - v) / 4.0,
//            (1 + u) * (1 - v) / 4.0,
//            (1 + u) * (1 + v) / 4.0,
//            (1 - u) * (1 + v) / 4.0,
//    };
    sf = {
            (1 - u) * (1 - v) / 4.0,
            (u) * (1 - v) / 4.0,
            (u) * (v) / 4.0,
            (1 - u) * (v) / 4.0,
    };
}

/*
 * Function uses the fact that hexahedra are always uniquely tied with the nodes of tetrahedra.
 * Therefore, knowing the barycentric coordinates inside a tetrahedron
 * makes it possible to use them to determine, in which part of the tetrahedron and therefore
 * also inside which hexahedron the point is located.
 */
int LinearQuadrangles::locate_cell(const Point3 &point, const int cell_guess) const {
    static constexpr int n_quads_per_tri = 3;

    int tri = quads->to_tri(abs(cell_guess));
    tri = lintri->locate_cell(point, tri);
    int sign = 1;
    if (tri < 0) sign = -1;
    tri = abs(tri);

    // calculate barycentric coordinates for a point
    array<double,3> bcc;
    lintri->get_shape_functions(bcc, point, tri);

    // all the ratios below are == 1, if the point is exactly on the boundary of two quadrangles
    // and 0 or inf, if point is on the edge of a triangle.
    // if ratio is > 1, the point is closer to 1st corresponding node, if < 1, it's closer to 2nd node
    const double b1b2 = bcc[0] / bcc[1];
    const double b1b3 = bcc[0] / bcc[2];
    const double b2b3 = bcc[1] / bcc[2];

    // point inside a quad connected to 1st triangular node ?
    if (b1b2 >= 1 && b1b3 >= 1)
        return sign * (n_quads_per_tri * tri + 0);
    // point inside a quad connected to 2nd triangular node ?
    if (b1b2 <= 1 && b2b3 >= 1)
        return sign * (n_quads_per_tri * tri + 1);
    // point inside a quad connected to 3rd triangular node ?
    if (b1b3 <= 1 && b2b3 <= 1)
        return sign * (n_quads_per_tri * tri + 2);

    return -1;
}

template class InterpolatorCells<3> ;
template class InterpolatorCells<6> ;
template class InterpolatorCells<4> ;
template class InterpolatorCells<10>;
template class InterpolatorCells<8>;

} /* namespace femocs */
