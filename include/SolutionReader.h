/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include <Interpolator.h>
#include "Primitives.h"
#include "Medium.h"

using namespace std;
namespace femocs {

/** Class to extract solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader conctructors */
    SolutionReader();
    SolutionReader(Interpolator* interpolator);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes
     * @param medium    atoms to be interpolated
     * @param r_cut     smoothing region cut off radius; 0 or less turns smoothing off
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data
     * @param srt      sort input atoms spatially
     */
    const void interpolate(const Medium &medium, double r_cut, int component=0, bool srt=true);

    /** Interpolate solution on points using the solution on tetrahedral mesh nodes
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data */
    const void interpolate(int n_points, double* x, double* y, double* z, double r_cut, int component=0, bool srt=true);

    /** Interpolate electric field on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void export_elfield(int n_points, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void export_potential(int n_points, double* phi, int* flag);

    /** Export calculated electic field distribution to HOLMOD */
    const void export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

private:
    Interpolator* interpolator;     ///< data needed for interpolation
    vector<Solution> interpolation; ///< interpolated data

    /** Function to clean the result from peaks
     * The cleaner makes the histogram for given component of electric field and applies smoothing
     * for the atoms, where field has diverging values.
     * For example, if histogram looks like [10 7 2 4 1 4 2 0 0 2], then the field on the two atoms that
     * made up the entries to last bin will be replaced by the average field around those two atoms. */
    const void clean(const int coordinate, const double r_cut);

    const void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);

    const Solution get_average_solution(const int I, const double r_cut);

    /** Reserve memory for interpolated data */
    const void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i) const;

    const void get_point_data(ofstream& out) const;
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
