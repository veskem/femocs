/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef LINEARINTERPOLATOR_H_
#define LINEARINTERPOLATOR_H_

#include "Primitives.h"
#include "Medium.h"
#include "TetgenMesh.h"
#include "Coarseners.h"
#include "laplace.h"
#include "currents_and_heating.h"

using namespace std;
namespace femocs {

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
 * */
class LinearInterpolator: public Medium {
public:
    /** LinearInterpolator conctructor */
    LinearInterpolator();

    /** Extract the electric potential and electric field values on the tetrahedra nodes from FEM solution */
    bool extract_solution(fch::Laplace<3>* laplace, const TetgenMesh &mesh);

    /** Extract the current density and temperature values on the tetrahedra nodes from FEM solution */
    bool extract_solution(fch::CurrentsAndHeatingStationary<3>* fem, const TetgenMesh &mesh);

    /** Interpolate both vector and scalar data.
     * Function assumes, that tetrahedron, that surrounds the point, is previously already found with locate_element.
     * @param point  point where the interpolation is performed
     * @param elem   tetrahedron around which the interpolation is performed */
    Solution get_solution(const Point3 &point, const int elem) const;

    /** Return the i-th entry from solution vector */
    Solution get_solution(const int i) const;

    /** Interpolate vector data.
     * Function assumes, that tetrahedron, that surrounds the point, is previously already found with locate_element.
     * @param point  point where the interpolation is performed
     * @param elem   tetrahedron around which the interpolation is performed */
    double get_scalar(const Point3 &point, const int elem) const;

    /** Return the i-th scalar entry from solution vector */
    double get_scalar(const int i) const;

    /** Interpolate scalar data.
     * Function assumes, that tetrahedron, that surrounds the point, is previously already found with locate_element.
     * @param point  point where the interpolation is performed
     * @param elem   tetrahedron around which the interpolation is performed */
    Vec3 get_vector(const Point3 &point, const int elem) const;

    /** Return the i-th Vec3 entry from solution vector */
    Vec3 get_vector(const int i) const;

    /** Locate the tetrahedron that surrounds or is closest to the point of interest.
     * The search starts from the elem_guess-th tetrahedron, proceedes with its neighbours
     * (number of neighbouring layers is specified inside the function),
     * then, if no match found, checks sequentially all the tetrahedra and if still no match found,
     * returns the index of tetrahedron whose centroid is closest to the point.
     * @param point       point of interest
     * @param elem_guess  index of tetrahedron around which the search starts
     * @return index of the tetrahedron that surrounds or is closest to the point
     */
    int locate_element(const Point3 &point, const int elem_guess);

    /** Print statistics about solution on node points */
    void print_statistics() const;

    /** Print the deviation from the analytical solution of hemiellipsoid on the infinite surface */
    void print_error(const Coarseners& c) const;

    /** Compare the analytical and calculated field enhancement */
    void print_enhancement() const;

    /** Set parameters to calculate analytical solution */
    void set_analyt(const Point3& origin, const double E0, const double radius1, const double radius2=-1);

    /** Enable or disable the search of points slightly outside the tetrahedron */
    void search_outside(const bool enable);

    /** Electric field that is assigned to atoms not found from mesh.
     *  Its value is BIG to make it immediately visible from the dataset. */
    const double error_field = 1e20;

private:
    /** Constants specifying the interpolation tolerances.
     * Making zero a bit negative allows searching points outside the tetrahedra. */
    const double epsilon = 0.1;
    double zero = -1.0 * epsilon;
    double radius1;  ///< Minor semi-axis of ellipse
    double radius2;  ///< Major semi-axis of ellipse
    double E0;       ///< Long-range electric field strength
    Point3 origin;

    vector<Solution> solution;          ///< interpolation data
    vector<SimpleElement> tetrahedra;   ///< tetrahedra node indices
    vector<vector<int>> tetneighbours;  ///< tetrahedra nearest neighbours
    vector<Point3> centroid;            ///< tetrahedra centroid coordinates

    vector<double> det0;                ///< major determinant for calculating bcc-s
    vector<Vec4> det1;                  ///< minor determinants for calculating 1st bcc
    vector<Vec4> det2;                  ///< minor determinants for calculating 2nd bcc
    vector<Vec4> det3;                  ///< minor determinants for calculating 3rd bcc
    vector<Vec4> det4;                  ///< minor determinants for calculating 4th bcc
    vector<bool> tet_not_valid;         ///< co-planarities of tetrahedra

    /** Pre-compute data about tetrahedra to make interpolation faster */
    void precompute_tetrahedra(const TetgenMesh &mesh);

    /** Return the mapping between tetrahedral and hexahedral meshes;
     * -1 indicates that mapping for corresponding object was not found */
    void get_maps(vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs,
            dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh, const double eps);

    /** Force the solution on tetrahedral nodes to be the weighed average
     * of the solutions on its Voronoi cell nodes */
    bool average_tetnodes(const TetgenMesh &mesh);

    /** Get barycentric coordinates for a point inside i-th tetrahedron */
    Vec4 get_bcc(const Point3 &point, const int i) const;

    /** Get whether the point is located inside the i-th tetrahedron */
    bool point_in_tetrahedron(const Point3 &point, const int i);

    /** Function to calculate determinant of 3x3 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2);

    /** Function to calculate determinant of 3x3 matrix */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);

    /** Function to calculate determinant of 4x4 matrix which's last column consists of ones */
    double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4);

    /** Function to calculate determinant of 4x4 matrix */
    double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;

    /** Get scalar and vector data associated with atoms */
    void get_cell_data(ofstream& outfile) const;

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Reserve memory for pre-compute data */
    void reserve_precompute(const int N);

    /** Return analytical potential value for i-th point near the hemisphere
     * @param radius  radius of the hemisphere
     * @param E0      long range electric field around the hemisphere */
    double get_analyt_potential(const int i, const Point3& origin) const;

    /** Return analytical electric field value for i-th point near the hemisphere
     * @param radius  radius of the hemisphere
     * @param E0      long range electric field around the hemisphere */
    Vec3 get_analyt_field(const int i) const;

    /** Get calculated field enhancement */
    double get_enhancement() const;

    /** Get analytical field enhancement for hemi-ellipsoid on infinite surface */
    double get_analyt_enhancement() const;
};

} // namespace femocs

#endif /* LINEARINTERPOLATOR_H_ */
