/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

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
class Interpolator: public Medium {
public:
    /** Interpolator conctructor */
    Interpolator();

    /** Extract the electric potential and electric field values on the tetrahedra nodes from FEM solution */
    const void extract_solution(DealII &fem, const TetgenMesh &mesh);

    const void extract_solution(fch::Laplace<3>* laplace, const TetgenMesh &mesh);
    
    const void extract_solution(fch::CurrentsAndHeating<3>* fem, const TetgenMesh &mesh);

    const Solution get_solution(const Point3 &point, const int elem);

    const double get_scalar(const Point3 &point, const int elem);

    const Vec3 get_vector(const Point3 &point, const int elem);

    const int locate_element(const Point3 &point, const int elem_guess);

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
    const void precompute_tetrahedra(const TetgenMesh &mesh);

    /** Return the mapping between tetrahedral and hexahedral meshes; -1 indicates that mapping for corresponding object was not found
     * @param fem         solution from Deal.II
     * @param tet2hex     mapping between tet- & hexmesh elements,
     * @param node2hex    mapping between tetmesh nodes & hexmesh elements,
     * @param node2vert   mapping between tetmesh nodes & hexmesh element's vertices. */
    const void get_maps(DealII& fem, vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert);

    const void get_maps(dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh,
            vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert);

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
    const string get_data_string(const int i) const;

    /** Reserve memory for interpolation data */
    const void reserve(const int N);

    /** Reserve memory for pre-compute data */
    const void reserve_precompute(const int N);

};

} // namespace femocs

#endif /* INTERPOLATOR_H_ */
