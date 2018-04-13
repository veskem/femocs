/*
 * InterpolatorCells.h
 *
 *  Created on: 10.1.2018
 *      Author: veske
 */

#ifndef INTERPOLATORCELLS_H_
#define INTERPOLATORCELLS_H_

#include "Globals.h"
#include "Primitives.h"
#include "TetgenMesh.h"
#include "TetgenCells.h"

using namespace std;
namespace femocs {

/**
 * Data & operations for obtaining and holding nodal interpolation data
 */
class InterpolatorNodes {
public:
    InterpolatorNodes();
    InterpolatorNodes(const string &norm_label, const string &scalar_label);
    ~InterpolatorNodes() {};

    /** Return number of available nodes */
    int size() const { return markers.size(); }

    /** Pre-compute data about vertices to make interpolation faster */
    void precompute();

    /** Output interpolation data in .xyz format */
    void write(const string &file_name) const;

    /** Output interpolation data to be appended to .vtk file */
    void write_point_data(ofstream& out) const;

    /** Print statistics about solution on node points */
    void print_statistics() const;

    /** Add solution to solutions vector */
    void append_solution(const Solution& solution) {
        expect(solutions.size() < solutions.capacity(), "Allocated vector size exceeded!");
        solutions.push_back(solution);
    }

    /** Return the pointer to the solution vector */
    vector<Solution>* get_solutions() { return &solutions; }

    /** Return full solution on i-th node */
    Solution get_solution(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return solutions[i];
    }

    /** Return the i-th vertex */
    Point3 get_vertex(const int i) const { return mesh->nodes[i]; }

    /** Return vector component of solution on i-th node */
    Vec3 get_vector(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return solutions[i].vector;
    }

    /** Return vector norm of solution on i-th node */
    double get_vector_norm(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return solutions[i].norm;
    }

    /** Return scalar component of solution on i-th node */
    double get_scalar(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return solutions[i].scalar;
    }

    /** Modify solution on the i-th node */
    void set_solution(const int i, const Solution& s) {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        solutions[i] = s;
    }

    /** Modify vector component of solution on the i-th node */
    void set_vector(const int i, const Vec3& v) {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        solutions[i].vector = v;
        solutions[i].norm = v.norm();
    }

    /** Modify scalar component of solution on the i-th node */
    void set_scalar(const int i, const double d) {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        solutions[i].scalar = d;
    }

    /** Change the data labels */
    void set_labels(const string& nl, const string& sl) {
        const_cast<string&>(norm_label) = nl;
        const_cast<string&>(scalar_label) = sl;
    }

    /** Change pointer to mesh */
    void set_mesh(const TetgenMesh* m) {
        mesh = m;
    }

    /** Return the max value of solution.norm data */
    double max_norm() const;

private:
    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc
    const string norm_label;        ///< description label attached to solution.norm -values
    const string scalar_label;      ///< description label attached to solution.scalar -values

    vector<Solution> solutions;     ///< interpolation data
    vector<int> markers;            ///< markers for nodes

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Output interpolation data in .xyz format */
    void write_xyz(ofstream& out) const;

    /** Output interpolation data in .vtk format */
    void write_vtk(ofstream& out) const;

    /** Return the cell type in vtk format */
    int get_cell_type() const { return TYPES.VTK.VERTEX; };
};

/**
 * Template class for interpolators of different kind
 */
template<int dim>
class InterpolatorCells {
public:
    InterpolatorCells() : mesh(NULL), nodes(NULL) { reserve(0); }

    InterpolatorCells(const InterpolatorNodes* n) :
        mesh(NULL), nodes(n) { reserve(0); }

    virtual ~InterpolatorCells() {};

    /** Return number of available cells */
    int size() const { return markers.size(); }

    /** Pick the suitable write function based on the file type.
     * Function is active only when file write is enabled */
    void write(const string &file_name) const;

    /** Pre-compute data about cells to make interpolation faster */
    virtual void precompute() {};

    /** Check whether the point is inside the cell */
    virtual bool point_in_cell(const Vec3& point, const int cell) const { return false; };

    /** Get interpolation weights for a point inside the cell */
    virtual array<double,dim> shape_functions(const Vec3& point, const int cell) const {
        require(false, "shape_functions(point, cell) not implemented for dim-" + d2s(dim));
        return array<double,dim>();
    }

    /** Get gradient of shape function for a point inside the cell */
    virtual array<Vec3,dim> shape_fun_grads(const Vec3& point, const int cell) const {
        require(false, "shape_fun_grads(point, cell) not implemented for dim-" + d2s(dim));
        return array<Vec3,dim>();
    }

    /** Get gradient of shape function for a cell node */
    virtual array<Vec3,dim> shape_fun_grads(const int cell, const int node) const {
        require(false, "shape_fun_grads(cell, node) not implemented for dim-" + d2s(dim));
        return array<Vec3,dim>();
    }

    /** Find the cell which contains the point or is the closest to it */
    virtual int locate_cell(const Point3 &point, const int cell_guess) const;

    /** @brief Interpolate both vector and scalar data inside or near the cell.
     * Function assumes that cell, that surrounds the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of shape functions and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed */
    virtual Solution interp_solution(const Point3 &point, const int cell) const;
    Solution interp_solution_v2(const Point3 &point, const int cell) const;

    /** Interpolate minus gradient of solution for any point inside a given cell
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed */
    Vec3 interp_gradient(const Point3 &point, const int cell) const;

    /** Interpolate minus gradient of solution for a node of a cell
     * @param cell   index of cell where the interpolation is performed
     * @param node   index of cell vertex where the interpolation is performed */
    Vec3 interp_gradient(const int cell, const int node) const;

    /** First locate the cell that surrounds or is closest to the point and then interpolate there.
     * Search starts from the argument cell. */
    Solution locate_interpolate(const Point3 &point, int& cell) const;
    Solution locate_interpolate_v2(const Point3 &point, int& cell) const;

    /** Modify cell marker */
    void set_marker(const int i, const int m) {
        require(i >= 0 && i < markers.size(), "Invalid index: " + d2s(i));
        markers[i] = m;
    }

    /** Access cell marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < markers.size(), "Invalid index: " + d2s(i));
        return markers[i];
    }

    /** Return i-th cell */
    virtual SimpleCell<dim> get_cell(const int i) const { return SimpleCell<dim>(); }

    /** Change the pointer to mesh */
    void set_mesh(const TetgenMesh* m) { mesh = m; }

    double decay_factor = -1.0;        ///< exp(decay_factor * node1.distance(node2)) gives the weight that can be used in smoothing process

protected:
    static constexpr double zero = 1e-15; ///< tolerance of calculations

    const TetgenMesh* mesh;           ///< Full mesh data with nodes, faces, elements etc
    const InterpolatorNodes* nodes;   ///< Pointer to nodes and solutions

    vector<int> markers;            ///< markers for cells
    vector<Point3> centroids;       ///< cell centroid coordinates
    vector<vector<int>> neighbours; ///< nearest neighbours of the cells

    /** Reserve memory for interpolation data */
    virtual void reserve(const int N);

    /** Return the cell type in vtk format */
    virtual int get_cell_type() const { return 0; };

    /** Output interpolation data in .vtk format */
    void write_vtk(ofstream& out) const;

    /** Output interpolation data to be appended to .vtk file */
    virtual void write_cell_data(ofstream& out) const;

    /** Find the common entry between two vectors */
    inline int common_entry(vector<unsigned>& vec1, vector<unsigned>& vec2) const;

    /** Determinant of 3x3 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2) const;

    /** Determinant of 3x3 matrix which's columns consist of Vec3-s */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const;

    /** Determinant of 4x4 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) const;

    /** Determinant of 4x4 matrix which's columns consist of Vec4-s */
    double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) const;
};

/**
 * Data & operations for linear tetrahedral interpolation without nodal data
 */
class LinearTetrahedra : public InterpolatorCells<4> {
public:
    LinearTetrahedra();
    LinearTetrahedra(const InterpolatorNodes* n);
    ~LinearTetrahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Get whether the point is located inside the i-th tetrahedron */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Get interpolation weights for a point inside i-th tetrahedron */
    array<double,4> shape_functions(const Vec3& point, const int i) const;

    /** Calculate the gradient of shape functions for a point inside i-th tetrahedron */
    array<Vec3, 4> shape_fun_grads(const Vec3& point, const int i) const;

    /** Return i-th tetrahedron */
    SimpleCell<4> get_cell(const int i) const { return (*tets)[i]; }

    /** Change the dependency data */
    void set_mesh(const TetgenMesh* m) {
        InterpolatorCells<4>::set_mesh(m);
        tets = &m->tets;
    }

    /** Specify the region where the cells are searched during the cell location. */
    void narrow_search_to(const int region);

private:
    const TetgenElements* tets;    ///< pointer to tetrahedra to access their specific routines

    vector<double> det0;            ///< major determinant for calculating bcc-s
    vector<Vec4> det1;              ///< minor determinants for calculating 1st bcc
    vector<Vec4> det2;              ///< minor determinants for calculating 2nd bcc
    vector<Vec4> det3;              ///< minor determinants for calculating 3rd bcc
    vector<Vec4> det4;              ///< minor determinants for calculating 4th bcc
    vector<bool> tet_not_valid;     ///< co-planarities of tetrahedra

    /** Reserve memory for pre-compute data */
    void reserve(const int N);

    /** Return the tetrahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.TETRAHEDRON; }
};

/**
 * Data & operations for quadratic tetrahedral interpolation without nodal data.
 * Class uses the data pre-computed for LinearTetrahedra.
 */
class QuadraticTetrahedra : public InterpolatorCells<10> {
public:
    QuadraticTetrahedra();
    QuadraticTetrahedra(const InterpolatorNodes* n, const LinearTetrahedra* l);
    ~QuadraticTetrahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3& point, const int cell) const;

    /** Find the tetrahedron which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Get interpolation weights for a point inside the tetrahedron */
    array<double,10> shape_functions(const Vec3& point, const int tet) const;
    
    /** Calculate the gradient of shape functions for a point inside the tetrahedron */
    array<Vec3,10> shape_fun_grads(const Vec3& point, const int tet) const;

    array<Vec3,10> shape_fun_grads_slow(const Vec3& point, const int tet) const;

    void test_shape_funs();

    /** Return i-th tetrahedron */
    SimpleCell<10> get_cell(const int i) const {
        require(i >= 0 && i < cells.size(), "Invalid index: " + d2s(i));
        return cells[i];
    }

    /** Change the dependency data */
    void set_mesh(const TetgenMesh* m) {
        InterpolatorCells<10>::set_mesh(m);
        tets = &m->tets;
    }

private:
    const TetgenElements* tets;    ///< pointer to tetrahedra to access their specific routines
    const LinearTetrahedra* lintet;   ///< Pointer to linear tetrahedra
    vector<QuadraticTet> cells;    ///< stored 10-noded tetrahedra

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the 10-noded tetrahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TETRAHEDRON; };

    /** Calculate the vertex indices of 10-noded tetrahedron */
    SimpleCell<10> calc_cell(const int i) const;
};

/**
 * Data & operations for linear hexahedral interpolation without nodal data.
 * Class uses the data pre-computed for LinearTetrahedra.
 */
class LinearHexahedra : public InterpolatorCells<8> {
public:
    LinearHexahedra();
    LinearHexahedra(const InterpolatorNodes* n, const LinearTetrahedra* l);
    ~LinearHexahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3 &point, const int cell) const;

    /** Find the hexahedron which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Get interpolation weights for a point inside i-th hexahedron */
    array<double,8> shape_functions(const Vec3& point, const int i) const;

    /** Get interpolation weights for a point inside i-th hexahedron
     *  and sort the result according to Deal.II ordering */
    array<double,8> shape_funs_dealii(const Vec3& point, const int i) const;

    /** Calculate the gradient of shape functions for a point inside i-th hexahedron */
    array<Vec3,8> shape_fun_grads(const Vec3& point, const int i) const;

    /** Calculate gradient of shape function for a hexahedral node */
    array<Vec3,8> shape_fun_grads(const int hex, const int node) const;

    /** Calculate the gradient of shape functions for a point inside the hexahedron
     * and sort the result according to Deal.II ordering */
    array<Vec3,8> shape_fun_grads_dealii(const Vec3& point, const int hex) const;

    /** Return the index of hexahedron in Deal.II that corresponds to i-th hexahedron;
     * -1 means there's no correspondence between two meshes */
    int femocs2deal(const int i) const {
        require(i >= 0 && i < map_femocs2deal.size(), "Invalid index: " + d2s(i));
        return map_femocs2deal[i];
    }

    /** Return the index of hexahedron in femocs that corresponds to i-th hexahedron in Deal.II */
    int deal2femocs(const int i) const {
        require(i >= 0 && i < map_deal2femocs.size(), "Invalid index: " + d2s(i) + ", size = " + d2s(map_deal2femocs.size()));
        return map_deal2femocs[i];
    }

    /** Return i-th hexahedron */
    SimpleCell<8> get_cell(const int i) const { return (*hexs)[i]; }

    /** Change the mesh dependency data */
    void set_mesh(const TetgenMesh* m) {
        InterpolatorCells<8>::set_mesh(m);
        hexs = &m->hexs;
    }

private:
    static constexpr int n_newton_iterations = 20; ///< max # Newton iterations while calculating natural coordinates

    const Hexahedra* hexs;          ///< pointer to hexahedra to access their specific routines
    const LinearTetrahedra* lintet; ///< Pointer to linear tetrahedra

    vector<int> map_femocs2deal;    ///< data for mapping between femocs and deal.II hex meshes
    vector<int> map_deal2femocs;    ///< data for mapping between deal.II and femocs hex meshes

    /// data for mapping point from Cartesian coordinates to natural ones
    vector<Vec3> f0s;
    vector<Vec3> f1s;
    vector<Vec3> f2s;
    vector<Vec3> f3s;
    vector<Vec3> f4s;
    vector<Vec3> f5s;
    vector<Vec3> f6s;
    vector<Vec3> f7s;

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return linear hexahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.HEXAHEDRON; };

    /** Map the point from Cartesian xyz-coordinates to natural uvw-coordinates.
     * In natural coordinate system, each coordinate is within limits [-1, 1]. */
    void project_to_nat_coords(double &u, double &v, double &w, const Vec3& point, const int hex) const;

    /** Look up for the natural coordinates u,v,w and nearest neighbouring nodes n1,n2,n3
     *  for the node inside a hexahedron. */
    void project_to_nat_coords(double &u, double &v, double &w,
        int &n1, int &n2, int &n3, const int node) const;
};

/**
 * Data & operations for linear triangular interpolation without nodal data
 * Class uses the data pre-computed for LinearTetrahedra.
 */
class LinearTriangles : public InterpolatorCells<3> {
public:
    LinearTriangles();
    LinearTriangles(const InterpolatorNodes* n, const LinearTetrahedra* lintet);
    ~LinearTriangles() {};

    /** Pre-compute data about triangles to make interpolation faster */
    void precompute();

    /** Check whether the projection of a point is inside the i-th triangle */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Calculate shape functions for a point with respect to the i-th triangle */
    array<double,3> shape_functions(const Vec3& point, const int i) const;

    /** Locate tetrahedron the surrounds the point and interpolate there.
     * Search starts from tetrahedra that are connected to the given triangle. */
    Solution interp_solution(const Point3 &point, const int tri) const;

    /** Interpolate conserved scalar data for the vector of atoms */
    void interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) const;

    /** Determine whether the point is within r_cut distance from the triangular surface.
     * If yes, the index of closest triangle will be returned. If no, -1 will be returned. */
    int near_surface(const Vec3& point, const double r_cut) const;

    /** Return the distance between a point and i-th triangle in the direction of its norm
     * without checking whether the projection of the point is inside the triangle or not */
    double fast_distance(const Vec3& point, const int i) const;

    /** Return norm of i-th triangle */
     Vec3 get_norm(const int i) const {
         require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
         return norms[i];
     }

     /** Return i-th triangle */
     SimpleCell<3> get_cell(const int i) const { return (*tris)[i]; }

     /** Change the mesh */
     void set_mesh(const TetgenMesh* m) {
         InterpolatorCells<3>::set_mesh(m);
         tris = &m->tris;
     }

private:
    const TetgenFaces* tris;         ///< Direct pointer to mesh triangles
    const LinearTetrahedra* lintet;  ///< Pointer to linear tetrahedron interpolator

    /** Data computed before starting looping through the triangles */
    vector<Vec3> vert0;
    vector<Vec3> edge1;
    vector<Vec3> edge2;
    vector<Vec3> pvec;
    vector<Vec3> norms;
    vector<double> max_distance;

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the triangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.TRIANGLE; };

    /** Return the distance between a point and i-th triangle in the direction of its norm.
     * If the projection of the point is outside the triangle, the 1e100 distance will be returned.
     * The calculations are based on Moller-Trumbore algorithm. The theory about it can be found from
     * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ */
    double distance(const Vec3& point, const int i) const;

    /** Append triangle specific data to vtk file */
    void write_cell_data(ofstream& out) const;
};

/**
 * Data & operations for quadratic triangular interpolation without nodal data.
 * Class uses the data pre-computed for LinearTriangles and QuadraticTetrahedra.
 */
class QuadraticTriangles : public InterpolatorCells<6> {
public:
    QuadraticTriangles();
    QuadraticTriangles(const InterpolatorNodes* n, const LinearTriangles* lintri, const QuadraticTetrahedra* quadtet);
    ~QuadraticTriangles() {};

    /** Pre-compute data about triangles to make interpolation faster */
    void precompute();

    /** Check whether the projection of a point is inside the i-th triangle */
    bool point_in_cell(const Vec3& point, const int cell) const;

    /** Find the triangle which contains the point projection or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Get interpolation weights for a point inside i-th triangle */
    array<double,6> shape_functions(const Vec3& point, const int i) const;

    /** Locate tetrahedron the surrounds the point and interpolate there.
     * Search starts from tetrahedra that are connected to the given triangle. */
    Solution interp_solution(const Point3 &point, const int tri) const;

    /** Return i-th hexahedron */
    SimpleCell<6> get_cell(const int i) const {
        require(i >= 0 && i < cells.size(), "Invalid index: " + d2s(i));
        return cells[i];
    }

    /** Change the mesh */
    void set_mesh(const TetgenMesh* m) {
        InterpolatorCells<6>::set_mesh(m);
        tris = &m->tris;
    }

private:
    const TetgenFaces* tris;            ///< Direct pointer to mesh triangles
    const LinearTriangles* lintri;      ///< Pointer to linear triangular interpolator
    const QuadraticTetrahedra* quadtet; ///< Pointer to quadratic tetrahedral interpolator
    vector<QuadraticTri> cells;         ///< stored 6-noded triangles

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the 6-noded triangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TRIANGLE; };

    /** Calculate the vertex indices of 6-noded triangle */
    SimpleCell<6> calc_cell(const int i) const;
};

/**
 * Data & operations for linear quadrangular interpolation without nodal data.
 * Class uses the data pre-computed for LinearTriangles and LinearHexahedra.
 */
class LinearQuadrangles : public InterpolatorCells<4> {
public:
    LinearQuadrangles();
    LinearQuadrangles(const InterpolatorNodes* n, const LinearTriangles* lintri, const LinearHexahedra* linhex);
    ~LinearQuadrangles() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3 &point, const int cell) const;

    /** Find the quadrangle which contains the point projection or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Locate hexahedron the surrounds the point and interpolate there.
     * Search starts from hexahedra that are connected to the given quadrangle. */
    Solution interp_solution(const Point3 &point, const int quad) const;

    /** Return i-th quadrangle */
    SimpleCell<4> get_cell(const int i) const { return (*quads)[i]; }

    /** Change the mesh */
    void set_mesh(const TetgenMesh* m) {
        InterpolatorCells<4>::set_mesh(m);
        quads = &m->quads;
    }

private:
    const Quadrangles* quads;       ///< Direct pointer to mesh quadrangles
    const LinearTriangles* lintri;  ///< Pointer to linear triangular interpolator
    const LinearHexahedra* linhex;  ///< Pointer to linear hexahedral interpolator

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the quadrangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRANGLE; };
};

} /* namespace femocs */

#endif /* INTERPOLATORCELLS_H_ */
