/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef LINEARINTERPOLATOR_H_
#define LINEARINTERPOLATOR_H_

#include "Primitives.h"
#include "TetgenMesh.h"
#include "TetgenCells.h"
#include "Coarseners.h"
#include "laplace.h"
#include "currents_and_heating.h"
#include "currents_and_heating_stationary.h"

using namespace std;
namespace femocs {

/** General class to linearly interpolate solution inside the mesh */
template<int dim>
class LinearInterpolator {

public:
    /** LinearInterpolator conctructor */
    LinearInterpolator() : vector_label("vector"), norm_label("vector_norm"), scalar_label("scalar"), mesh(NULL), nodes(NULL) {
        reserve(0);
        reserve_precompute(0);
    }

    LinearInterpolator(const TetgenMesh* m) : vector_label("vector"), norm_label("vector_norm"),
            scalar_label("scalar"), mesh(m), nodes(&m->nodes) {
        reserve(0);
        reserve_precompute(0);
    }

    LinearInterpolator(const TetgenMesh* m, const string& vl, const string& nl, const string& sl) :
        vector_label(vl), norm_label(nl), scalar_label(sl), mesh(m), nodes(&m->nodes) {
        reserve(0);
        reserve_precompute(0);
    }

    virtual ~LinearInterpolator() {};

    /** Return number of available interpolation nodes */
    int size() const { return solutions.size(); }

    /** Pick the suitable write function based on the file type.
     * Function is active only when file write is enabled */
    void write(const string &file_name) const;

    /** Enable or disable the search of points slightly outside the cell */
    void search_outside(const bool enable) {
        if (enable) {
            zero = -1.0 * epsilon;
            one = 1.0 + epsilon;
        } else {
            zero = 0;
            one = 1.0;
        }
    }

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

    /** Return scalar component of solution on i-th node */
    double get_scalar(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return solutions[i].scalar;
    }

    /** @brief Interpolate both vector and scalar data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed; >= 0 for cells around the point and < 0 for others */
    Solution interp_solution(const Point3 &point, const int cell) const;

    /** @brief Interpolate vector data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed; >= 0 for cells around the point and < 0 for others */
    Vec3 interp_vector(const Point3 &point, const int cell) const;

    /** @brief Interpolate scalar data inside or near the cell.
     * Function assumes, that cell, that fits the best to the point, is previously already found with locate_cell.
     * cell>=0 initiates the usage of barycentric coordinates and cell<0 the usage of mere distance-dependent weighting.
     * @param point  point where the interpolation is performed
     * @param cell   index of cell around which the interpolation is performed; >= 0 for cells around the point and < 0 for others */
    double interp_scalar(const Point3 &point, const int cell) const;

    /** Find the cell which contains the point or is the closest to it */
    int locate_cell(const Point3 &point, const int cell_guess);

protected:
    /** Constants specifying the interpolation tolerances.
     * Making zero a bit negative allows searching points outside the tetrahedra. */
    const double epsilon = 0.1;
    double zero = -1.0 * epsilon;
    double one = 1.0 + epsilon;
    double decay_factor = -1.0;     ///< exp(decay_factor * node1.distance(node2)) gives the weight that can be used in smoothing process

    const string vector_label;      ///< description label attached to solution.vector -values
    const string norm_label;        ///< description label attached to solution.norm -values
    const string scalar_label;      ///< description label attached to solution.scalar -values

    vector<Solution> solutions;     ///< interpolation data
    vector<vector<int>> neighbours; ///< nearest neighbours of the cells
    vector<Point3> centroids;       ///< cell centroid coordinates
    vector<Point3> vertices;        ///< coordinates of cell vertices
    vector<SimpleCell<dim>> cells;  ///< interpolation cells - triangles or tetrahedra

    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc
    const TetgenNodes* nodes;       ///< Mesh nodes

    /** Reserve memory for interpolation data */
    void reserve(const int N) {
        require(N >= 0, "Invalid number of points: " + to_string(N));
        solutions.clear();
        solutions.reserve(N);
    }

    /** Reserve memory for pre-computation data */
    virtual void reserve_precompute(const int N) {
        neighbours = vector<vector<int>>(N);
        centroids.clear();
        centroids.reserve(N);
        cells.clear();
        cells.reserve(N);
        vertices.clear();
        vertices.reserve(nodes->size());
    }

    /** Pre-compute data about cells to make interpolation faster */
    virtual void precompute() {}

    /** Return the cell type in vtk format;
     * 5-triangle, 9-quadrangle, 10-tetrahedron, 12-hexahedron */
    virtual int get_cell_type() const { return 0; }

    /** Calculate barycentric coordinates for a point with respect to the cell */
    virtual array<double,dim> get_bcc(const Vec3& point, const int cell) const { return array<double,dim>(); }

    /** Calculate distance-dependent weights for a point with respect to the cell */
    array<double,dim> get_weights(const Point3 &point, const SimpleCell<dim>& scell) const;

    /** Check whether the point is inside the cell.
     * It does not use get_bcc routine to achieve faster performance. */
    virtual bool point_in_cell(const Vec3& point, const int cell) { return false; }

    /** Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
     *  surrounding hexahedral nodes */
    bool average_sharp_nodes(const vector<vector<unsigned>>& vorocells);

    /** Return the mapping between tetrahedral and hexahedral meshes;
     * -1 indicates that mapping for corresponding object was not found */
    void get_maps(vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs,
            dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh, const double eps);

    /** Output node data in .xyz format */
    void write_xyz(ofstream& out, const int n_nodes) const;

    /** Output interpolation cell data in .vtk format */
    void write_vtk(ofstream& out, const int n_nodes) const;

    /** Determinant of 3x3 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2) {
        return v1.x * (v2.y - v2.z) - v1.y * (v2.x - v2.z) + v1.z * (v2.x - v2.y);
    }

    /** Determinant of 3x3 matrix which's columns consist of Vec3-s */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) {
        return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z)
                + v3.x * (v1.y * v2.z - v2.y * v1.z);
    }

    /** Determinant of 4x4 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4) {
        const double det1 = determinant(v2, v3, v4);
        const double det2 = determinant(v1, v3, v4);
        const double det3 = determinant(v1, v2, v4);
        const double det4 = determinant(v1, v2, v3);

        return det4 - det3 + det2 - det1;
    }

    /** Determinant of 4x4 matrix which's columns consist of Vec4-s */
    double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4) {
        double det1 = determinant(Vec3(v1.y,v1.z,v1.w), Vec3(v2.y,v2.z,v2.w), Vec3(v3.y,v3.z,v3.w));
        double det2 = determinant(Vec3(v1.x,v1.z,v1.w), Vec3(v2.x,v2.z,v2.w), Vec3(v3.x,v3.z,v3.w));
        double det3 = determinant(Vec3(v1.x,v1.y,v1.w), Vec3(v2.x,v2.y,v2.w), Vec3(v3.x,v3.y,v3.w));
        double det4 = determinant(Vec3(v1.x,v1.y,v1.z), Vec3(v2.x,v2.y,v2.z), Vec3(v3.x,v3.y,v3.z));

        return v4.w * det4 - v4.z * det3 + v4.y * det2 - v4.x * det1;
    }
};

/** Class to linearly interpolate solution inside tetrahedral mesh */
/* Useful links
 *
 * Compact theory how to find barycentric coordintes (bbc):
 * http://steve.hollasch.net/cgindex/geometry/ptintet.html
 *
 * Properties of determinant:
 * http://www.vitutor.com/alg/determinants/properties_determinants.html
 * http://www.vitutor.com/alg/determinants/minor_cofactor.html
 *
 * c++ code to find and handle bcc-s:
 * http://dennis2society.de/painless-tetrahedral-barycentric-mapping
 *
 * Interpolating inside the element using bcc:
 * http://www.cwscholz.net/projects/diss/html/node37.html
 *
 */
class TetrahedronInterpolator : public LinearInterpolator<4> {
public:
    /** TetrahedronInterpolator conctructor */
    TetrahedronInterpolator(const TetgenMesh* mesh);

    /** Extract the electric potential and electric field values on the tetrahedra nodes from FEM solution */
    bool extract_solution(fch::Laplace<3>* laplace);

    /** Extract the current density and temperature values on the tetrahedra nodes from FEM solution */
    bool extract_solution(fch::CurrentsAndHeatingStationary<3>* fem);

    /** Extract the current density and temperature values on the tetrahedra nodes from FEM solution */
    bool extract_solution(fch::CurrentsAndHeating<3>* fem);

    /** Print statistics about solution on node points */
    void print_statistics() const;

private:
    const TetgenElements* elems;    ///< Direct pointer to tetrahedra to access their specific routines

    vector<double> det0;            ///< major determinant for calculating bcc-s
    vector<Vec4> det1;              ///< minor determinants for calculating 1st bcc
    vector<Vec4> det2;              ///< minor determinants for calculating 2nd bcc
    vector<Vec4> det3;              ///< minor determinants for calculating 3rd bcc
    vector<Vec4> det4;              ///< minor determinants for calculating 4th bcc
    vector<bool> tet_not_valid;     ///< co-planarities of tetrahedra

    /** Force the solution on tetrahedral nodes to be the weighed average
     * of the solutions on its Voronoi cell nodes */
    bool average_sharp_nodes(const bool vacuum);

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute();

    /** Reserve memory for pre-compute data */
    void reserve_precompute(const int N);

    /** Get whether the point is located inside the i-th tetrahedron */
    bool point_in_cell(const Vec3& point, const int i);

    /** Get barycentric coordinates for a point inside i-th tetrahedron */
    array<double,4> get_bcc(const Vec3& point, const int i) const;

    /** Return the tetrahedron type in vtk format */
    int get_cell_type() const { return 10; }
};

/** Class to interpolate solution on surface triangles */
class TriangleInterpolator : public LinearInterpolator<3> {
public:
    /** Constructor of TriangleInterpolator  */
    TriangleInterpolator(const TetgenMesh* mesh);

    /** Extract the electric potential and electric field values on triangular mesh nodes from FEM solution */
    bool extract_solution(fch::Laplace<3>* fem);

    /** Calculate charges on surface faces using direct solution in the face centroids */
    void calc_charges(const double E0);

    /** Interpolate conserved scalar data for the vector of atoms */
    void interp_conserved(vector<double>& scalars, const vector<Atom>& atoms);

    /** Print statistics about solution on node points */
    void print_statistics(const double Q);

    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]

private:
    const TetgenFaces* faces;    ///< Direct pointer to triangles to access their specific routines

    /** Data computed before starting looping through the triangles */
    vector<Vec3> vert0;
    vector<Vec3> edge1;
    vector<Vec3> edge2;
    vector<Vec3> pvec;
    vector<Vec3> norms;
    vector<double> max_distance;

    bool clean_nodes();

    /** Force the solution on triangular nodes to be the weighed average
     * of the solutions on its surrounding quadrangular nodes */
    bool average_sharp_nodes(const bool vacuum);

    /** Precompute the data needed to calculate the distance of points from surface
     * in the direction of triangle surface norms */
    void precompute();

    /** Reserve memory for precompute data */
    void reserve_precompute(const int n);

    /** Check whether the projection of a point is inside the i-th triangle */
    bool point_in_cell(const Vec3& point, const int i);

    /** Calculate barycentric coordinates for a point with respect to the i-th triangle */
    array<double,3> get_bcc(const Vec3& point, const int i) const;

    /** Return the triangle type in vtk format */
    int get_cell_type() const { return 5; }
};

} // namespace femocs

#endif /* LINEARINTERPOLATOR_H_ */
