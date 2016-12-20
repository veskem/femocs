/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include "Primitives.h"
#include "Medium.h"
#include "Interpolator.h"

using namespace std;
namespace femocs {

/** Class to extract solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader conctructors */
    SolutionReader();
    SolutionReader(Interpolator* interpolator);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void interpolate(const Medium &medium);

    /** Interpolate electric field on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void export_elfield(int n_points, double* x, double* y, double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    const void export_potential(int n_points, double* x, double* y, double* z, double* phi, int* flag);

    /** Export calculated electic field distribution to HOLMOD */
    const void export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Function to clean the result from peaks
     * The cleaner makes the histogram for given component of electric field and applies smoothing
     * for the atoms, where field has diverging values.
     * For example, if histogram looks like [10 7 2 4 1 4 2 0 0 2], then the field on the two atoms that
     * made up the entries to last bin will be replaced by the average field around those two atoms. */
    const void clean(const int coordinate, const int n_bins, const double smooth_factor, const double r_cut);

    /** Output atom data in .vtk format */
    const void write_vtk(const string &file_name) const;

private:
    Interpolator* interpolator;     ///< data needed for interpolation
    vector<Solution> interpolation; ///< interpolated data
    
    const void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);

    const Solution get_average_solution(const int I, const double smooth_factor, const double r_cut);

    /** Reserve memory for interpolated data */
    const void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i) const;
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
