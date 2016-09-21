/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include "Primitives.h"
#include "DealII.h"

using namespace std;
namespace femocs {

/** Class to extract solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader conctructor */
    SolutionReader();

    /** Extract the electric potential and electric field values on Medium atoms from FEM solution */
    const void extract_solution(DealII* fem, Medium &medium);
    const void extract_statistics(Mesh &mesh);
    const void smoothen_result(const double smooth_width);
    const void export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    const void print_statistics();
    const void sort_atoms(const int x1, const int x2, const string& direction = "up");

    const Solution get_solution(const int i);

private:
    DealII* fem;
    double longrange_efield;
    vector<double> face_qualities;
    vector<double> elem_qualities;
    vector<Solution> solution;

    /** Electric field that is assigned to atoms not found from mesh.
     *  Its value is BIG to make it immediately visible from data set. */
    const double error_field = 1e20;

    const void smoothen_result_ema(const double smooth_width);
    const void smoothen_result_sma(const int n_samples);

    inline Vec3 get_ema(const int i0, const int i1, const double smooth_width);
    inline Vec3 get_sma_up(const int i, const int n_samples);
    inline Vec3 get_sma_down(const int i, const int n_samples);

    /** Reserve memory for solution vectors */
    const void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);

    const vector<int> get_node2face_map(Mesh &mesh, int node);
    const vector<int> get_node2elem_map(Mesh &mesh, int node);

    /**
     * Return the mapping between atoms & nodes, nodes & elements and nodes & vertices.
     * In medium2node the value -1 indicates that there's no node in the mesh that corresponds to the given atom.
     */
    const void get_maps(Medium& medium, vector<int>& medium2node, vector<int>& node2elem, vector<int>& node2vert);
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
 * */
class Interpolator {
public:
    /** Interpolator conctructor */
    Interpolator();
    Interpolator(Mesh* mesh, SolutionReader* solution);

    const void reserve(const int N);

    const void test();

private:
    SolutionReader* solution;
    Mesh* mesh;

    vector<Point3> centroid;
    vector<double> det0;
    vector<Vec4> det1;
    vector<Vec4> det2;
    vector<Vec4> det3;
    vector<Vec4> det4;

    const void precompute_tetrahedra();

    const double determinant(const Vec3 &v1, const Vec3 &v2);
    const double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);
    const double determinant(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, const Vec3 &v4);
    const double determinant(const Vec4 &v1, const Vec4 &v2, const Vec4 &v3, const Vec4 &v4);
    const Vec4 get_bcc(const Point3 &point, const int elem);
    const Vec3 get_interpolation(const Point3 &point, const int elem);

};

} /* namespace femocs */

#endif /* SOLUTIONREADER_H_ */
