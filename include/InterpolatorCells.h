/*
 * InterpolatorCells.h
 *
 *  Created on: 10.1.2018
 *      Author: veske
 */

#ifndef INTERPOLATORCELLS_H_
#define INTERPOLATORCELLS_H_

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
    InterpolatorNodes(const TetgenMesh* m, const string& norm_label, const string& scalar_label);
    ~InterpolatorNodes() {};

    /** Return number of available nodes */
    int size() const { return vertices.size(); }

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
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i];
    }

    /** Modify solution on the i-th node */
    void set_solution(const int i, const Solution& s) {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        solutions[i] = s;
    }

    /** Return the pointer to the vertices vector */
    vector<Point3>* get_vertices() { return &vertices; }

    /** Return the i-th vertex */
    Point3 get_vertex(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return vertices[i];
    }

    /** Return vector component of solution on i-th node */
    Vec3 get_vector(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i].vector;
    }

    /** Return vector norm of solution on i-th node */
    double get_vector_norm(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i].norm;
    }

    /** Return scalar component of solution on i-th node */
    double get_scalar(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i].scalar;
    }

    /** Change the dependency data */
    void set_dependencies(const TetgenMesh* m, const string& nl, const string& sl) {
        mesh = const_cast<TetgenMesh*>(m);
        const_cast<string&>(norm_label) = nl;
        const_cast<string&>(scalar_label) = sl;
    }

private:
    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc
    const string norm_label;        ///< description label attached to solution.norm -values
    const string scalar_label;      ///< description label attached to solution.scalar -values

    vector<Solution> solutions;     ///< interpolation data
    vector<Point3> vertices;        ///< coordinates of cell vertices
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

    InterpolatorCells(const TetgenMesh* m, const InterpolatorNodes* n) :
        mesh(m), nodes(n) { reserve(0); }

    virtual ~InterpolatorCells() {};

    /** Return number of available cells */
    int size() const { return cells.size(); }

    /** Pick the suitable write function based on the file type.
     * Function is active only when file write is enabled */
    void write(const string &file_name) const;

    /** Pre-compute data about cells to make interpolation faster */
    virtual void precompute() {};

    /** Check whether the point is inside the cell */
    virtual bool point_in_cell(const Vec3& point, const int cell) const { return false; };

    /** Get interpolation weights for a point inside i-th tetrahedron */
    virtual void get_shape_functions(array<double,dim>& sf, const Vec3& point, const int i) const {}

    /** Find the cell which contains the point or is the closest to it */
    virtual int locate_cell(const Point3 &point, const int cell_guess) const;

    /** @brief Interpolate both vector and scalar data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed */
    Solution interp_solution(const Point3 &point, const int c) const;

    /** Change the dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n) {
        mesh = const_cast<TetgenMesh*>(m);
        nodes = const_cast<InterpolatorNodes*>(n);
    }

    /** Modify cell marker */
    void set_marker(const int i, const int m) {
        require(i >= 0 && i < markers.size(), "Invalid index: " + to_string(i));
        markers[i] = m;
    }

    /** Access cell marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < markers.size(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    double decay_factor = -1.0;        ///< exp(decay_factor * node1.distance(node2)) gives the weight that can be used in smoothing process

protected:
    static constexpr double zero = 1e-15; ///< tolerance of calculations

    const TetgenMesh* mesh;           ///< Full mesh data with nodes, faces, elements etc
    const InterpolatorNodes* nodes;   ///< Pointer to nodes and solutions

    vector<int> markers;            ///< markers for cells
    vector<Point3> centroids;       ///< cell centroid coordinates
    vector<vector<int>> neighbours; ///< nearest neighbours of the cells
    vector<SimpleCell<dim>> cells;  ///< interpolation cells

    /** Reserve memory for interpolation data */
    virtual void reserve(const int N);

    /** Return the cell type in vtk format */
    virtual int get_cell_type() const { return 0; };

    /** Calculate distance-dependent weights for a point with respect to the cell */
    void get_weights(array<double,dim>& weights, const Point3 &point, const SimpleCell<dim>& scell) const;

    /** Output interpolation data in .vtk format */
    void write_vtk(ofstream& out) const;

    /** Output interpolation data to be appended to .vtk file */
    virtual void write_cell_data(ofstream& out) const;

    /** Find the common entry between two vectors */
    inline int common_entry(vector<unsigned>& vec1, vector<unsigned>& vec2) const;
};

/**
 * Data & operations for linear tetrahedral interpolation without nodal data
 */
class LinearTetrahedra : public InterpolatorCells<4> {
public:
    LinearTetrahedra();
    LinearTetrahedra(const TetgenMesh* m, const InterpolatorNodes* n);
    ~LinearTetrahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Get whether the point is located inside the i-th tetrahedron */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Get interpolation weights for a point inside i-th tetrahedron */
    void get_shape_functions(array<double,4>& sf, const Vec3& point, const int i) const;

    /** Change the dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n) {
        InterpolatorCells<4>::set_dependencies(m, n);
        elems = &m->elems;
    }

    /** Specify the region where the cells are searched during the cell location. */
    void narrow_search_to(const int region);

    /** Determinant of 3x3 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2) const;

    /** Determinant of 3x3 matrix which's columns consist of Vec3-s */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const;

    /** Determinant of 4x4 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) const;

    /** Determinant of 4x4 matrix which's columns consist of Vec4-s */
    double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) const;

private:
    const TetgenElements* elems;    ///< pointer to tetrahedra to access their specific routines

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
    QuadraticTetrahedra(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTetrahedra* l);
    ~QuadraticTetrahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3& point, const int cell) const;

    /** Get interpolation weights for a point inside i-th tetrahedron */
    void get_shape_functions(array<double,10>& sf, const Vec3& point, const int i) const;

    /** Change the dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTetrahedra* l) {
        InterpolatorCells<10>::set_dependencies(m, n);
        elems = &m->elems;
        lintet = const_cast<LinearTetrahedra*>(l);
    }

private:
    const TetgenElements* elems;    ///< pointer to tetrahedra to access their specific routines
    const LinearTetrahedra* lintet;   ///< Pointer to linear tetrahedra

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the 10-noded tetrahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TETRAHEDRON; };

    /** Calculate the vertex indices of 10-noded tetrahedron */
    SimpleCell<10> get_cell(const int i) const;
};

/**
 * Data & operations for linear hexahedral interpolation without nodal data.
 * Class uses the data pre-computed for LinearTetrahedra.
 */
class LinearHexahedra : public InterpolatorCells<8> {
public:
    LinearHexahedra();
    LinearHexahedra(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTetrahedra* l);
    ~LinearHexahedra() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Get interpolation weights for a point inside i-th hexahedron */
    void get_shape_functions(array<double,8>& sf, const Vec3& point, const int i) const;

    /** Find the hexahedron which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Return the index of hexahedron in Deal.II that corresponds to i-th hexahedron;
     * -1 means there's no correspondence between two meshes */
    int femocs2deal(const int i) const {
        require(i >= 0 && i < map_femocs2deal.size(), "Invalid index: " + to_string(i));
        return map_femocs2deal[i];
    }

    /** Return the index of hexahedron in femocs that corresponds to i-th hexahedron in Deal.II */
    int deal2femocs(const int i) const {
        require(i >= 0 && i < map_deal2femocs.size(), "Invalid index: " + to_string(i) + "size = " + to_string(map_deal2femocs.size()));
        return map_deal2femocs[i];
    }

    /** Change the mesh dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTetrahedra* l) {
        InterpolatorCells<8>::set_dependencies(m, n);
        hexs = &m->hexahedra;
        lintet = const_cast<LinearTetrahedra*>(l);
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

    /** Return the linear hexahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.HEXAHEDRON; };

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3& point, const int cell) const {
        require(false, "point_in_cell not used in LinearHexahedra!");
        return false;
    }
};

/**
 * Data & operations for linear triangular interpolation without nodal data
 */
class LinearTriangles : public InterpolatorCells<3> {
public:
    LinearTriangles();
    LinearTriangles(const TetgenMesh* m, const InterpolatorNodes* n);
    ~LinearTriangles() {};

    /** Pre-compute data about triangles to make interpolation faster */
    void precompute();

    /** Check whether the projection of a point is inside the i-th triangle */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Calculate shape functions for a point with respect to the i-th triangle */
    void get_shape_functions(array<double,3>& sf, const Vec3& point, const int i) const;

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
         require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
         return norms[i];
     }

     /** Change the dependency data */
     void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n) {
         InterpolatorCells<3>::set_dependencies(m, n);
         faces = &m->faces;
     }

private:
    const TetgenFaces* faces;    ///< Direct pointer to triangles to access their specific routines

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
 * Class uses the data pre-computed for LinearTriangles.
 */
class QuadraticTriangles : public InterpolatorCells<6> {
public:
    QuadraticTriangles();
    QuadraticTriangles(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTriangles* l);
    ~QuadraticTriangles() {};

    /** Pre-compute data about triangles to make interpolation faster */
    void precompute();

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3& point, const int cell) const;

    /** Get interpolation weights for a point inside i-th triangle */
    void get_shape_functions(array<double,6>& sf, const Vec3& point, const int i) const;

    /** Change the dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTriangles* l) {
        InterpolatorCells<6>::set_dependencies(m, n);
        faces = &m->faces;
        lintri = const_cast<LinearTriangles*>(l);
    }

private:
    const TetgenFaces* faces;    ///< Direct pointer to triangles to access their specific routines
    const LinearTriangles* lintri;   ///< Pointer to linear triangles

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the 6-noded triangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TRIANGLE; };

    /** Calculate the vertex indices of 6-noded triangle */
    SimpleCell<6> get_cell(const int i) const;
};

/**
 * Data & operations for linear quadrangular interpolation without nodal data.
 * Class uses the data pre-computed for LinearTriangles.
 */
class LinearQuadrangles : public InterpolatorCells<4> {
public:
    LinearQuadrangles();
    LinearQuadrangles(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTriangles* l);
    ~LinearQuadrangles() {};

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Get interpolation weights for a point inside i-th hexahedron */
    void get_shape_functions(array<double,4>& sf, const Vec3& point, const int i) const;

    /** Find the hexahedron which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Until proper way to calculate shape functions is found, use lintri interpolator */
    Solution interp_solution(const Point3 &point, const int c) const;

    /** Change the mesh dependency data */
    void set_dependencies(const TetgenMesh* m, const InterpolatorNodes* n, const LinearTriangles* l) {
        InterpolatorCells<4>::set_dependencies(m, n);
        quads = &m->quads;
        lintri = const_cast<LinearTriangles*>(l);
    }

private:
    const Quadrangles* quads;       ///< pointer to quadrangles to access their specific routines
    const LinearTriangles* lintri;  ///< Pointer to linear triangles

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Return the quadrangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRANGLE; };

    /** Check whether the point is inside the cell */
    bool point_in_cell(const Vec3& point, const int cell) const {
        require(false, "point_in_cell not used in LinearQuadrangles!");
        return false;
    }
};

} /* namespace femocs */

#endif /* INTERPOLATORCELLS_H_ */
