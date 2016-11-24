/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include "Primitives.h"
#include "SolutionReader.h"
#include "Medium.h"
#include "TetgenMesh.h"

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
class Interpolator: public Medium {
public:
    /** Interpolator conctructor */
    Interpolator();
    Interpolator(SolutionReader* sr);

    /** Pre-compute data about tetrahedra to make interpolation faster */
    const void precompute_tetrahedra(const TetgenMesh &mesh);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void extract_interpolation(const Medium &medium);

    /** Interpolate electric field on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void extract_elfield(int n_points, double* x, double* y, double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void extract_potential(int n_points, double* x, double* y, double* z, double* phi, int* flag);

    /** Export calculated electic field distribution to HOLMOD */
    const void export_interpolation(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Function to clean the result from peaks
     * The cleaner makes the histogram for given component of electric field and applies smoothing
     * for the atoms, where field has diverging values.
     * For example, if histogram looks like [10 7 2 4 1 4 2 0 0 2], then the field on the two atoms that
     * made up the entries to last bin will be replaced by the average field around those two atoms. */
    const void clean(const int coordinate, const int n_bins, const double smooth_factor, const double r_cut);

private:
    /** Constants specifying the interpolation tolerances.
     * Making zero a bit negative allows to interpolate outside the tetrahedron. */
    const double epsilon = 1e-1;
    const double zero = -1.0 * epsilon;

    SolutionReader* solution;           ///< solution data
    vector<Solution> interpolation;     ///< interpolation data

    vector<SimpleElement> tetrahedra;   ///< tetrahedra node indices
    vector<vector<int>> tetneighbours;  ///< tetrahedra nearest neighbours
    vector<Point3> centroid;            ///< tetrahedra centroid coordinates

    vector<double> det0;                ///< major determinant for calculating bcc-s
    vector<Vec4> det1;                  ///< minor determinants for calculating 1st bcc
    vector<Vec4> det2;                  ///< minor determinants for calculating 2nd bcc
    vector<Vec4> det3;                  ///< minor determinants for calculating 3rd bcc
    vector<Vec4> det4;                  ///< minor determinants for calculating 4th bcc
    vector<bool> tet_not_valid;         ///< co-planarities of tetrahedra

    const Solution get_interpolation(const Point3 &point, const int elem);
    const void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);
    const Solution get_average_solution(const int I, const double smooth_factor, const double r_cut);
    const int locate_element(const Point3 &point, const int elem_guess);
    const Vec4 get_bcc(const Point3 &point, const int elem);
    const bool point_in_tetrahedron(const Point3 &point, const int i);

    /** Function to calculate determinant of 3x3 matrix which's last column consists of ones */
    const double determinant(const Vec3 &v1, const Vec3 &v2);

    /** Function to calculate determinant of 3x3 matrix */
    const double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);

    /** Function to calculate determinant of 4x4 matrix which's last column consists of ones */
    const double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4);

    /** Function to calculate determinant of 4x4 matrix */
    const double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);

    /** Reserve memory for interpolation data */
    const void reserve(const int N);

    /** Reserve memory for pre-compute data */
    const void reserve_precompute(const int N);
};

} // namespace femocs

#endif /* INTERPOLATOR_H_ */
