/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include "LinearInterpolator.h"
#include "Primitives.h"
#include "Medium.h"
#include "TetgenCells.h"
#include "TetgenMesh.h"

using namespace std;
namespace femocs {

/** Class to interpolate solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader constructors */
    SolutionReader();
    SolutionReader(LinearInterpolator* ip, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab);

protected:
    const string vec_label;       ///< label for vector data
    const string vec_norm_label;  ///< label for data associated with vector length
    const string scalar_label;    ///< label for scalar data

    LinearInterpolator* interpolator;     ///< data needed for interpolation
    vector<Solution> interpolation; ///< interpolated data

    /** Reserve memory for data */
    void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors for .xyz file;
     * i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;

    /** Get information about data vectors for .vtk file. */
    void get_point_data(ofstream& out) const;
};

/** Class to extract solution from DealII calculations */
class FieldReader: public SolutionReader {
public:
    /** SolutionReader constructors */
    FieldReader();
    FieldReader(LinearInterpolator* ip);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes
     * @param medium    atoms to be interpolated
     * @param r_cut     smoothing region cut off radius; 0 or less turns smoothing off
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data
     * @param srt      sort input atoms spatially
     */
    void interpolate(const Medium &medium, const double r_cut, const int component=0, const bool srt=true);

    /** Interpolate solution on points using the solution on tetrahedral mesh nodes
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data */
    void interpolate(const int n_points, const double* x, const double* y, const double* z,
            const double r_cut, const int component=0, const bool srt=true);

    /** Interpolate electric field on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    void export_elfield(const int n_points, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    void export_potential(const int n_points, double* phi, int* flag);

    /** Export calculated electic field distribution to HOLMOD */
    void export_solution(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Print statistics about interpolated solution */
    void print_statistics();

private:
    /** Function to clean the result from peaks
     * The cleaner makes the histogram for given component of electric field and applies smoothing
     * for the atoms, where field has diverging values.
     * For example, if histogram looks like [10 7 2 4 1 4 2 0 0 2], then the field on the two atoms that
     * made up the entries to last bin will be replaced by the average field around those two atoms. */
    void clean(const int coordinate, const double r_cut);

    void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);

    Solution get_average_solution(const int I, const double r_cut);
};

/** Class to interpolate current densities and temperatures */
class HeatReader: public SolutionReader {
public:
    /** HeatReader constructors */
    HeatReader();
    HeatReader(LinearInterpolator* ip);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes */
    void interpolate(const Medium &medium, const double r_cut);

private:
};

/** Class to calculate charges from electric field */
class ChargeReader: public SolutionReader {
public:
    /** ChargeReader constructors */
    ChargeReader();
    ChargeReader(LinearInterpolator* ip);

    void calc_interpolated_charges(const TetgenMesh& mesh);

    void calc_charges(const TetgenMesh& mesh);

    Vec3 get_force(const int i) const;

    double get_area(const int i) const;

    double get_charge(const int i) const;

    void print_statistics(const Medium::Sizes& sizes, const double E0);

private:
    const double eps0 = 8.854187817620e-12; ///< vacuum permittivity [F/m]

    void get_elfields(const TetgenMesh& mesh, vector<Vec3> &elfields);

    double get_total_charge() const;
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
