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
#include "DealII.h"
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
    void extract_solution(fch::Laplace<3>* laplace, const TetgenMesh &mesh, const double w0_elfield, const double w0_pot);
    
    void extract_solution(DealII* fem, const TetgenMesh &mesh, const double w0_elfield, const double w0_pot);

    void extract_solution(fch::CurrentsAndHeating<3>* fem, const TetgenMesh &mesh);

    Solution get_solution(const Point3 &point, const int elem);

    double get_scalar(const Point3 &point, const int elem);

    Vec3 get_vector(const Point3 &point, const int elem);

    int locate_element(const Point3 &point, const int elem_guess);

    /** Electric field that is assigned to atoms not found from mesh.
     *  Its value is BIG to make it immediately visible from the dataset. */
    const double error_field = 1e20;

private:
    /** Constants specifying the interpolation tolerances.
     * Making zero a bit negative allows to interpolate outside the tetrahedron. */
    const double epsilon = 1e-1;
    const double zero = -1.0 * epsilon;

    vector<Solution> solution;     ///< interpolation data

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
    void get_maps(dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh,
            vector<int>& tet2hex, vector<int>& cell_indxs, vector<int>& vert_indxs);

    /** Force the solution on tetrahedral nodes to be the weighed average
     * of the solutions on its Voronoi cell nodes */
    void average_tetnodes(const TetgenMesh &mesh, const double w0_elfield, const double w0_pot);

    /** Get barycentric coordinates for a point inside i-th tetrahedron */
    Vec4 get_bcc(const Point3 &point, const int i);

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

    /** Reserve memory for interpolation data */
    void reserve(const int N);

    /** Reserve memory for pre-compute data */
    void reserve_precompute(const int N);

};

} // namespace femocs

#endif /* LINEARINTERPOLATOR_H_ */
