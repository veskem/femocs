/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include "Primitives.h"
#include "TetgenMesh.h"
#include "TetgenCells.h"
#include "InterpolatorCells.h"
#include "Coarseners.h"
#include "laplace.h"
#include "currents_and_heating.h"
#include "currents_and_heating_stationary.h"

using namespace std;
namespace femocs {

/** General class for interpolating solution inside mesh
 * Class holds subclasses (InterpolatorCells) that hold the actual interpolation data and routines
 * that allow performing interpolation either on surface or in space.
 * 
 * For tetrahedral and triangular interpolation, barycentric coordinates (BCC-s)
 * are used to calculate the shape (aka interpolation) functions.
 * Compact theory how to find BCC-s:
 * http://steve.hollasch.net/cgindex/geometry/ptintet.html
 *
 * Properties of determinant:
 *   http://www.vitutor.com/alg/determinants/properties_determinants.html
 *   http://www.vitutor.com/alg/determinants/minor_cofactor.html
 *
 * c++ code to find and handle BCC-s:
 *   http://dennis2society.de/painless-tetrahedral-barycentric-mapping
 *
 * Theory about shape functions inside different primitives.
 * Linear & quadratic triangle & quadrangle:
 *   https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch18.d/IFEM.Ch18.pdf
 * Linear tetrahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
 * Quadratic tetrahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
 * Linear hexahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
 */
class Interpolator {
public:
    Interpolator(const TetgenMesh* m, const string& norm_label, const string& scalar_label);
    ~Interpolator() {};

    /** Extract the electric potential and field values from FEM solution */
    bool extract_solution(fch::Laplace<3>* fem);

    /** Extract the current density and stationary temperature values from FEM solution */
    bool extract_solution(fch::CurrentsAndHeatingStationary<3>* fem);

    /** Extract the current density and transient temperature values from FEM solution */
    bool extract_solution(fch::CurrentsAndHeating<3>& fem);

    InterpolatorNodes nodes;      ///< vertices and solutions on them
    LinearTriangles lintris;      ///< data & operations for linear triangular interpolation
    QuadraticTriangles quadtris;  ///< data & operations for quadratic triangular interpolation
    LinearTetrahedra lintets;     ///< data & operations for linear tetrahedral interpolation
    QuadraticTetrahedra quadtets; ///< data & operations for quadratic tetrahedral interpolation
    LinearHexahedra linhexs;      ///< data & operations for linear hexahedral interpolation

private:
    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc

    /** Calculate the mapping between Femocs & deal.II mesh nodes,
     *  nodes & hexahedral elements and nodes & element's vertices.
     *  -1 indicates that mapping for corresponding object was not found */
    void get_maps(vector<int>& femocs2deal, vector<int>& cell_indxs, vector<int>& vert_indxs,
            dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh) const;

    /** Transfer solution from FEM solver to Interpolator */
    void store_solution(const vector<int>& femocs2deal,
            const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data);

    /** Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
     *  surrounding hexahedral nodes */
    bool average_sharp_nodes(const bool vacuum);
};

/** General class for interpolating solution inside mesh
 * To make it possible to choose between different types of interpolators,
 * the following is inspired by the Strategy design pattern.
 *
 * Useful links about Strategy design pattern:
 * https://en.wikipedia.org/wiki/Strategy_pattern
 * https://r3dux.org/2011/07/an-example-implementation-of-the-strategy-design-pattern-in-c/
 */
class GeneralInterpolator {
public:
    GeneralInterpolator() :
        mesh(NULL), nodes(NULL), norm_label("vector_norm"), scalar_label("scalar")
    {
        reserve(0);
        reserve_precompute(0);
    }

    GeneralInterpolator(const TetgenMesh* m) :
        mesh(m), nodes(&m->nodes), norm_label("vector_norm"), scalar_label("scalar")
    {
        reserve(0);
        reserve_precompute(0);
    }

    GeneralInterpolator(const TetgenMesh* m, const string& nl, const string& sl) :
        mesh(m), nodes(&m->nodes), norm_label(nl), scalar_label(sl)
    {
        reserve(0);
        reserve_precompute(0);
    }

    virtual ~GeneralInterpolator() {};

    /** Return number of available interpolation nodes */
    int size() const { return solutions.size(); }

    /** Pick the suitable write function based on the file type.
     * Function is active only when file write is enabled */
    void write(const string &file_name) const;

    /** Add solution to solutions vector */
    void append_solution(const Solution& solution) {
        expect((unsigned)size() < solutions.capacity(), "Allocated vector size exceeded!");
        solutions.push_back(solution);
    }
    
    /** Return the pointer to the solution vector */
    vector<Solution>* get_solutions() { return &solutions; }

    /** Return full solution on i-th node */
    Solution get_solution(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i];
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

    /** Find the cell which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess) const;

    /** Find the cell which contains the point or is the closest to it looking only for non-marked cells */
    int locate_cell(const Point3 &point, const int cell_guess, vector<bool>& cell_checked) const;

    /** Extract the electric potential and field values from FEM solution */
    bool extract_solution(fch::Laplace<3>* fem);

    /** Extract the current density and stationary temperature values from FEM solution */
    bool extract_solution(fch::CurrentsAndHeatingStationary<3>* fem);

    /** Extract the current density and transient temperature values from FEM solution */
    bool extract_solution(fch::CurrentsAndHeating<3>* fem);

    /** Pre-compute data about cells to make interpolation faster */
    virtual void precompute() = 0;

    /** Interpolate both vector and scalar data inside or near the cell. */
    virtual Solution interp_solution(const Point3 &point, const int cell) const = 0;

protected:
    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]
    const double coplanar_epsilon = 0.1; ///< coplanarity tolerance
    double vertex_epsilon = 0;      ///< max distance between two identical vertices
    double decay_factor = -1.0;     ///< exp(decay_factor * node1.distance(node2)) gives the weight that can be used in smoothing process
    
    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc
    const TetgenNodes* nodes;       ///< Mesh nodes

    const string norm_label;        ///< description label attached to solution.norm -values
    const string scalar_label;      ///< description label attached to solution.scalar -values

    vector<Solution> solutions;     ///< interpolation data
    vector<Point3> centroids;       ///< cell centroid coordinates
    vector<Point3> vertices;        ///< coordinates of cell vertices
    vector<vector<int>> neighbours; ///< nearest neighbours of the cells
   
    /** Reserve memory for interpolation data */
    void reserve(const int N) {
        require(N >= 0, "Invalid number of points: " + to_string(N));
        solutions.clear();
        solutions.reserve(N);
    }

    /** Reserve memory for pre-computation data */
    virtual void reserve_precompute(const int N) {
        require(N >= 0, "Invalid number of points: " + to_string(N));
        require(nodes, "No mesh found!");

        centroids.clear();
        centroids.reserve(N);
        neighbours = vector<vector<int>>(N);
        vertices.clear();
        vertices.reserve(nodes->size());
    }

    /** Return the cell type in vtk format */
    virtual int get_cell_type() const = 0;

    /** Check whether the point is inside the cell */
    virtual bool point_in_cell(const Vec3& point, const int cell) const = 0;

    /** Output the indices of cell vertices */
    virtual void write_cells(ofstream& out) const = 0;

    /** Output interpolation cell data in .vtk format */
    virtual void write_vtk(ofstream& out) const;

    /** Output interpolation data in .xyz format */
    void write_xyz(ofstream& out) const;

    /** Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
     *  surrounding hexahedral nodes */
    virtual bool average_sharp_nodes(const bool vacuum) = 0;

    bool average_sharp_nodes(const vector<vector<unsigned>>& vorocells);

    /** Calculate the mapping between Femocs & deal.II mesh nodes,
     *  nodes & hexahedral elements and nodes & element's vertices.
     *  -1 indicates that mapping for corresponding object was not found */
    void get_maps(vector<int>& femocs2deal, vector<int>& cell_indxs, vector<int>& vert_indxs,
            dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh) const;

    /** Transfer solution from FEM solver to TemplateInterpolator */
    void store_solution(const vector<int>& femocs2deal,
            const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data);

    /** Find the common entry between two vectors */
    int common_entry(vector<unsigned>& vec1, vector<unsigned>& vec2) const {
        for (unsigned i : vec1)
            for (unsigned j : vec2)
                if (i == j) return i;
        return -1;
    }
};

/** General class for interpolating on the surface */
class SurfaceInterpolator : public GeneralInterpolator {
 public:
    SurfaceInterpolator(const TetgenMesh* m) : GeneralInterpolator(m) {}
    virtual ~SurfaceInterpolator() {};

    /** Function to determine the distance of a point from the surface */
    virtual double fast_distance(const Vec3& point, const int i) const = 0;

    /** Return norm of i-th triangle */
    virtual Vec3 get_norm(const int i) const = 0;

    /** Interpolate conserved scalar data for the vector of atoms */
    virtual void interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) const = 0;
};

/** General class for interpolating in the space */
class VolumeInterpolator : public GeneralInterpolator {
public:
    VolumeInterpolator(const TetgenMesh* m) : GeneralInterpolator(m) {}
    virtual ~VolumeInterpolator() {};
};

/** Class to interpolate solution inside tetrahedral mesh */
template<int rank>
class TetrahedronInterpolator : public VolumeInterpolator {
public:
    TetrahedronInterpolator(const TetgenMesh* m) : VolumeInterpolator(m), elems(&m->elems) {}
    virtual ~TetrahedronInterpolator() {};

    /** @brief Interpolate both vector and scalar data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed */
    Solution interp_solution(const Point3 &point, const int cell) const;

    /** Print statistics about solution on node points */
    void print_statistics() const;

protected:
    const TetgenElements* elems;    ///< pointer to tetrahedra to access their specific routines

    vector<double> det0;            ///< major determinant for calculating bcc-s
    vector<Vec4> det1;              ///< minor determinants for calculating 1st bcc
    vector<Vec4> det2;              ///< minor determinants for calculating 2nd bcc
    vector<Vec4> det3;              ///< minor determinants for calculating 3rd bcc
    vector<Vec4> det4;              ///< minor determinants for calculating 4th bcc
    vector<bool> tet_not_valid;     ///< co-planarities of tetrahedra

    vector<SimpleCell<rank>> cells;  ///< interpolation tetrahedra

    /** Force the solution on tetrahedral nodes to be the weighed average
     * of the solutions on its Voronoi cell nodes */
    bool average_sharp_nodes(const bool vacuum);

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Reserve memory for pre-compute data */
    void reserve_precompute(const int N);

    /** Get whether the point is located inside the i-th tetrahedron */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Get interpolation weights for a point inside i-th tetrahedron */
    virtual void get_shape_functions(array<double,rank>& sf, const Vec3& point, const int i) const {}

    /** Calculate distance-dependent weights for a point with respect to the cell */
    void get_weights(array<double,rank>& weights, const Point3 &point, const SimpleCell<rank>& scell) const;

    /** Return the tetrahedron type in vtk format */
    virtual int get_cell_type() const { return 0; }

    /** Get stored tetrahedron */
    virtual SimpleCell<rank> get_cell(const int i) const { return SimpleCell<rank>(); }
    
    /** Write vertex indices */
    void write_cells(ofstream& out) const {
        const int n_cells = cells.size();
        out << "\nCELLS " << n_cells << " " << (1+rank) * n_cells << "\n";
        for (int i = 0; i < n_cells; ++i)
            out << rank << " " << cells[i] << "\n";
    }

        /** Determinant of 3x3 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2) const {
        return v1.x * (v2.y - v2.z) - v1.y * (v2.x - v2.z) + v1.z * (v2.x - v2.y);
    }

    /** Determinant of 3x3 matrix which's columns consist of Vec3-s */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const {
        return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z)
                + v3.x * (v1.y * v2.z - v2.y * v1.z);
    }

    /** Determinant of 4x4 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) const {
        const double det1 = determinant(v2, v3, v4);
        const double det2 = determinant(v1, v3, v4);
        const double det3 = determinant(v1, v2, v4);
        const double det4 = determinant(v1, v2, v3);

        return det4 - det3 + det2 - det1;
    }

    /** Determinant of 4x4 matrix which's columns consist of Vec4-s */
    double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) const {
        double det1 = determinant(Vec3(v1.y,v1.z,v1.w), Vec3(v2.y,v2.z,v2.w), Vec3(v3.y,v3.z,v3.w));
        double det2 = determinant(Vec3(v1.x,v1.z,v1.w), Vec3(v2.x,v2.z,v2.w), Vec3(v3.x,v3.z,v3.w));
        double det3 = determinant(Vec3(v1.x,v1.y,v1.w), Vec3(v2.x,v2.y,v2.w), Vec3(v3.x,v3.y,v3.w));
        double det4 = determinant(Vec3(v1.x,v1.y,v1.z), Vec3(v2.x,v2.y,v2.z), Vec3(v3.x,v3.y,v3.z));

        return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
    }
};

/** Class to interpolate solution inside 4-noded tetrahedra */
class LinTetInterpolator : public TetrahedronInterpolator<4> {
public:
    LinTetInterpolator(const TetgenMesh* mesh) : TetrahedronInterpolator<4>(mesh) {}

private:
    /** Get interpolation weights for a point inside i-th tetrahedron */
     void get_shape_functions(array<double,4>& sf, const Vec3& point, const int i) const;

     /** Return the tetrahedron type in vtk format */
     int get_cell_type() const { return TYPES.VTK.TETRAHEDRON; }

     SimpleCell<4> get_cell(const int i) const;
};

/** Class to interpolate solution inside 10-noded tetrahedra */
class QuadTetInterpolator : public TetrahedronInterpolator<10> {
public:
    QuadTetInterpolator(const TetgenMesh* mesh) : TetrahedronInterpolator<10>(mesh) {}

private:
    /** Get interpolation weights for a point inside i-th tetrahedron */
    void get_shape_functions(array<double,10>& sf, const Vec3& point, const int i) const;

    /** Return the tetrahedron type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TETRAHEDRON; }
     
    SimpleCell<10> get_cell(const int i) const;

    /** Locate the index of node that is in the centroid of opposite face of the given tetrahedral vertex */
    int opposite_node(const int tet, const int vert) const;

    double integrate(const int hex) const;

    void hex_shape_function(array<double,8>& sf, const double u, const double v, const double w) const;

    bool average_and_check_sharp_nodes(const bool vacuum);
};

/** Class to interpolate solution inside surface triangles */
template<int rank>
class TriangleInterpolator : public SurfaceInterpolator {
public:
    TriangleInterpolator(const TetgenMesh* m) : SurfaceInterpolator(m), faces(&m->faces) {}
    virtual ~TriangleInterpolator() {};

    /** Print statistics about solution on node points */
    void print_statistics(const double Q);

    /** @brief Interpolate both vector and scalar data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed */
    Solution interp_solution(const Point3 &point, const int cell) const;

    /** Interpolate conserved scalar data for the vector of atoms */
    void interp_conserved(vector<double>& scalars, const vector<Atom>& atoms) const;

    /** Precompute the data needed to calculate the distance of points from surface
     * in the direction of triangle surface norms */
    void precompute();

    /** Determine whether the point is within r_cut distance from the triangular surface */
    int near_surface(const Vec3& point, const double r_cut) const;

    /** Return norm of i-th triangle */
    Vec3 get_norm(const int i) const {
        require(i >= 0 && i < this->size(), "Invalid index: " + to_string(i));
        return norms[i];
    }

    /** Return the distance between a point and i-th triangle in the direction of its norm */
    double fast_distance(const Vec3& point, const int i) const;

protected:
    const TetgenFaces* faces;    ///< Direct pointer to triangles to access their specific routines

    /** Data computed before starting looping through the triangles */
    vector<Vec3> vert0;
    vector<Vec3> edge1;
    vector<Vec3> edge2;
    vector<Vec3> pvec;
    vector<Vec3> norms;
    vector<double> max_distance;

    vector<SimpleCell<rank>> cells;  ///< interpolation triangles

    /** Add triangle specific data to vtk file */
    void write_vtk(ofstream& out) const;

    /** Write vertex indices */
    void write_cells(ofstream& out) const {
        const int n_cells = cells.size();
        out << "\nCELLS " << n_cells << " " << (1+rank) * n_cells << "\n";
        for (int i = 0; i < n_cells; ++i)
            out << rank << " " << cells[i] << "\n";
    }

    /** Force the solution on triangular nodes to be the weighed average
     * of the solutions on its surrounding quadrangular nodes */
    bool average_sharp_nodes(const bool vacuum);

    /** Reserve memory for precompute data */
    void reserve_precompute(const int n);

    /** Return the distance between a point and i-th triangle in the direction of its norm.
     * The calculations are based on Moller-Trumbore algorithm. The theory about it can be found from
     * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ */
    double distance(const Vec3& point, const int i) const;

    /** Check whether the projection of a point is inside the i-th triangle */
    bool point_in_cell(const Vec3& point, const int i) const;

    /** Calculate shape functions for a point with respect to the i-th triangle */
    virtual void get_shape_functions(array<double,rank>& sf, const Vec3& point, const int i) const {}

    /** Calculate distance-dependent weights for a point with respect to the cell */
    void get_weights(array<double,rank>& weights, const Point3 &point, const SimpleCell<rank>& scell) const;

    /** Return the triangle type in vtk format */
    virtual int get_cell_type() const { return 0; }

    /** Get stored triangle */
    virtual SimpleCell<rank> get_cell(const int i) const { return SimpleCell<rank>(); }
};

/** Class to interpolate solution inside 3-noded triangles */
class LinTriInterpolator : public TriangleInterpolator<3> {
public:
    LinTriInterpolator(const TetgenMesh* mesh) : TriangleInterpolator<3>(mesh) {}

private:
    /** Calculate shape functions for a point with respect to the i-th triangle */
    void get_shape_functions(array<double,3>& sf, const Vec3& point, const int i) const;

    /** Return the triangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.TRIANGLE; }

    /** Get stored triangle */
    SimpleCell<3> get_cell(const int i) const;
};

/** Class to interpolate solution inside 6-noded triangles */
class QuadTriInterpolator : public TriangleInterpolator<6> {
public:
    QuadTriInterpolator(const TetgenMesh* mesh) : TriangleInterpolator<6>(mesh) {}

private:
    /** Calculate shape functions for a point with respect to the i-th triangle */
    void get_shape_functions(array<double,6>& sf, const Vec3& point, const int i) const;

    /** Return the 6-noded triangle type in vtk format */
    int get_cell_type() const { return TYPES.VTK.QUADRATIC_TRIANGLE; }

    /** Get stored triangle */
    SimpleCell<6> get_cell(const int i) const;
};

} // namespace femocs

#endif /* INTERPOLATOR_H_ */
