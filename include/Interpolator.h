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
#include "Mesh.h"

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

    /** Interpolate in medium atoms using the solution on tetrahedral mesh nodes */
    const void extract_interpolation(SolutionReader* solution, const Medium &medium);

    /** Export calculated electic field distribution to HOLMOD */
    const void export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Function to clean the result from peaks */
    const void clean(const int coordinate, const int n_bins, const double smooth_factor, const double r_cut);

private:
    /** Constants to specify the tolerances.
     * Making zero a bit negative allows to interpolate outside the tetrahedron */
    const double epsilon = 1e-2;
    const double zero = -1.0 * epsilon;

    /** Electric field that is assigned to atoms not found from mesh.
     *  Its value is BIG to make it immediately visible from data set. */
    const double error_field = 1e20;

    Mesh* mesh;                      //!< tetrahedral mesh
    SolutionReader* solution;        //!< solution data
    vector<Solution> interpolation;  //!< interpolation data

    /** Vectors holding precomputed data */
    vector<Point3> centroid;
    vector<double> det0;
    vector<Vec4> det1;
    vector<Vec4> det2;
    vector<Vec4> det3;
    vector<Vec4> det4;
    vector<bool> tet_not_valid;

    const Solution get_interpolation(const Point3 &point, const int elem);
    const void precompute_tetrahedra();

    const void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);
    const Vec3 get_average_solution(const int I, const double smooth_factor, const double r_cut);
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

    /** Reserve memory for interpolation and pre-compute vectors */
    const void reserve(const int N);

};

} /* namespace femocs */

#endif /* INTERPOLATOR_H_ */
