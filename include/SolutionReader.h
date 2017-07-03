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
#include "currents_and_heating.h"
#include "getelec.h"

using namespace std;
namespace femocs {

/** Class to interpolate solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader constructors */
    SolutionReader();
    SolutionReader(LinearInterpolator* ip, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab);

    /** Reserve memory for data */
    void reserve(const int n_nodes);

    /** Print statistics about interpolated solution */
    void print_statistics();

    /** Compare interpolated scalar statistics with a constant */
    void print_statistics(const double Q);

    /** Append solution */
    void append_interpolation(const Solution& s);

    /** Get i-th Solution */
    Solution get_interpolation(const int i) const;

    /** Set i-th Solution */
    void set_interpolation(const int i, const Solution& s);

    /** Calculate statistics about coordinates and solution */
    void calc_statistics();

    /** Statistics about solution */
    struct Statistics {
        double vec_norm_min;  ///< minimum value of vector norm
        double vec_norm_max;  ///< maximum value of vector norm
        double scal_min;      ///< minimum value of scalar
        double scal_max;      ///< maximum value of scalar
    } stat;

protected:
    const string vec_label;       ///< label for vector data
    const string vec_norm_label;  ///< label for data associated with vector length
    const string scalar_label;    ///< label for scalar data

    LinearInterpolator* interpolator;  ///< data needed for interpolation
    vector<Solution> interpolation;    ///< interpolated data

    /** Initialise statistics about coordinates and solution */
    void init_statistics();

    /** Get i-th entry from all data vectors for .xyz file;
     * i < 0 gives the header of data vectors */
    string get_data_string(const int i) const;

    /** Get information about data vectors for .vtk file. */
    void get_point_data(ofstream& out) const;

    /** Function to clean the result from peaks
     * The cleaner makes the histogram for given component of electric field and applies smoothing
     * for the atoms, where field has diverging values.
     * For example, if histogram looks like [10 7 2 4 1 4 2 0 0 2], then the field on the two atoms that
     * made up the entries to last bin will be replaced by the average field around those two atoms. */
    void clean(const int coordinate, const double r_cut);

    void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);

    Solution get_average_solution(const int I, const double r_cut);
};

/** Class to extract solution from DealII calculations */
class FieldReader: public SolutionReader {
public:
    /** SolutionReader constructors */
    FieldReader();
    FieldReader(LinearInterpolator* ip);

    void interpolate(const double r_cut, const int component, const bool srt);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes
     * @param medium    atoms to be interpolated
     * @param r_cut     smoothing region cut off radius; 0 or less turns smoothing off
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data
     * @param srt       sort input atoms spatially
     */
    void interpolate(const Medium &medium, const double r_cut, const int component=0, const bool srt=true);

    /** Interpolate solution on points using the solution on tetrahedral mesh nodes
     * @param component component of result to interpolate: 0-all, 1-vector data, 2-scalar data */
    void interpolate(const int n_points, const double* x, const double* y, const double* z,
            const double r_cut, const int component=0, const bool srt=true);

    /** Calculate the electric field for the current and temperature solver */
    void interpolate(fch::CurrentsAndHeating<3>* ch_solver, const double r_cut, const bool srt=true);

    /** Interpolate electric field on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    void export_elfield(const int n_points, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    void export_potential(const int n_points, double* phi, int* flag);

    /** Export calculated electic field distribution to HOLMOD */
    void export_solution(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    Vec3 get_elfield(const int i) const;

    double get_potential(const int i) const;

    /** Compare the analytical and calculated field enhancement */
    void print_enhancement() const;

    /** Set parameters to calculate analytical solution */
    void set_analyt(const double E0, const double radius1, const double radius2=-1);

    void calc_emission(fch::CurrentsAndHeating<3>* ch_solver);
    
    void emission_line(const Point3& point, const Vec3& field,
                        vector<double> &rline, vector<double> &Vline, double rmax);

private:
    double radius1;  ///< Minor semi-axis of ellipse
    double radius2;  ///< Major semi-axis of ellipse
    double E0;       ///< Long-range electric field strength

    /** Get calculated field enhancement */
    double get_enhancement() const;

    /** Get analytical field enhancement for hemi-ellipsoid on infinite surface */
    double get_analyt_enhancement() const;
};

/** Class to interpolate current densities and temperatures */
class HeatReader: public SolutionReader {
public:
    /** HeatReader constructors */
    HeatReader();
    HeatReader(LinearInterpolator* ip);

    /** Interpolate solution on medium atoms using the solution on tetrahedral mesh nodes */
    void interpolate(const Medium &medium);

    /** Export interpolated temperature */
    void export_temperature(const int n_atoms, double* T);

    Vec3 get_rho(const int i) const;

    double get_temperature(const int i) const;

private:
};

/** Class to calculate charges from electric field */
class ChargeReader: public SolutionReader {
public:
    /** ChargeReader constructors */
    ChargeReader();
    ChargeReader(LinearInterpolator* ip);

    void calc_interpolated_charges(const TetgenMesh& mesh, const double E0);

    void calc_charges(const TetgenMesh& mesh, const double E0);

    void clean(const Medium::Sizes& sizes, const double latconst);

    Vec3 get_elfield(const int i) const;

    double get_area(const int i) const;

    double get_charge(const int i) const;

    /** Check whether charge is conserved within specified limits */
    bool charge_conserved(const double Q, const double eps) const;

    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]

private:
    void get_elfields(const TetgenMesh& mesh, vector<Vec3> &elfields);
};

/** Class to calculate forces from charges and electric fields */
class ForceReader: public SolutionReader {
public:
    /** ChargeReader constructors */
    ForceReader();
    ForceReader(LinearInterpolator* ip);

    /** Replace the charge and force on the nanotip nodes with the one found with Voronoi cells */
    void recalc_forces(const FieldReader &fields, const vector<Vec3>& areas, const double force_factor);
        
    /** Calculate forces from atomic electric fields and face charges */
    void calc_forces(const FieldReader &fields, const ChargeReader& faces,
        const double r_cut, const double smooth_factor, const double force_factor);

    /** Export the induced charge and force on imported atoms
     * @param n_atoms  number of first atoms field is calculated
     * @param xq       charge and force in PARCAS format (xq[0] = q1, xq[1] = Fx1, xq[2] = Fy1, xq[3] = Fz1, xq[4] = q2, xq[5] = Fx2 etc)
     */
    void export_force(const int n_atoms, double* xq);

    Vec3 get_force(const int i) const;

    double get_charge(const int i) const;

    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]

private:
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
