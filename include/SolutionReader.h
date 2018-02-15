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
#include "TetgenCells.h"
#include "TetgenMesh.h"
#include "VoronoiMesh.h"
#include "Interpolator.h"
#include "Config.h"
#include "CurrentsAndHeating.h"
#include "CurrentsAndHeatingStationary.h"

using namespace std;
namespace femocs {

/** Class to interpolate solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader constructors */
    SolutionReader();
    SolutionReader(Interpolator* i, const string& vec_lab, const string& vec_norm_lab, const string& scal_lab);

    /** Interpolate solution on the system atoms */
    void calc_interpolation();

    /** Interpolate solution using already available data about which atom is connected with which cell */
    void calc_interpolation(vector<int>& atom2face);

    /** Reserve memory for data */
    void reserve(const int n_nodes);

    /** Print statistics about interpolated solution */
    void print_statistics();

    /** Append solution */
    void append_interpolation(const Solution& s);

    /** Get pointer to interpolation vector */
    vector<Solution>* get_interpolations();

    /** Get i-th Solution */
    Solution get_interpolation(const int i) const;

    /** Set i-th Solution */
    void set_interpolation(const int i, const Solution& s);

    /** Set interpolation preferences */
    void set_preferences(const bool _srt, const int _dim, const int _rank) {
        require((_dim == 2 || _dim == 3), "Invalid interpolation dimension: " + to_string(_dim));
        require((_rank == 1 || _rank == 2 || _rank == 3), "Invalid interpolation rank: " + to_string(_rank));
        sort_atoms = _srt;
        dim = _dim;
        rank = _rank;
    }

    /** Alter the pointer to interpolator */
    void set_interpolator(Interpolator* i) { interpolator = i; }

    /** Calculate statistics about coordinates and solution */
    void calc_statistics();

    /** Check and replace the NaNs and peaks in the solution */
    bool clean(const double r_cut, const bool use_hist_clean);

    /** Statistics about solution */
    struct Statistics {
        double vec_norm_min;  ///< minimum value of vector norm
        double vec_norm_max;  ///< maximum value of vector norm
        double scal_min;      ///< minimum value of scalar
        double scal_max;      ///< maximum value of scalar
    } stat;

    int get_femocs_index(int deal_index){
        return interpolator->linhexs.deal2femocs(deal_index);
    }

    int export_results(const int n_points, const string &data_type, const bool append, double* data);

    int interpolate_results(const int n_points, const string &data_type, const double* x,
            const double* y, const double* z, double* data);

protected:
    const string vec_label;       ///< label for vector data
    const string vec_norm_label;  ///< label for data associated with vector length
    const string scalar_label;    ///< label for scalar data
    double limit_min;             ///< minimum value of accepted comparison value
    double limit_max;             ///< maximum value of accepted comparison value
    bool sort_atoms;              ///< sort atoms along Hilbert curve to make interpolation faster
    int dim;                      ///< location of interpolation; 2-surface, 3-space
    int rank;                     ///< interpolation rank; 1-linear, 2-quadratic

    Interpolator* interpolator;    ///< pointer to interpolator
    vector<Solution> interpolation;       ///< interpolated data

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
    void histogram_clean(const int coordinate, const double r_cut);

    void get_histogram(vector<int> &bins, vector<double> &bounds, const int coordinate);

    Solution get_average_solution(const int I, const double r_cut);
};

/** Class to extract solution from DealII calculations */
class FieldReader: public SolutionReader {
public:

    FieldReader(Interpolator* i);

    void test_pic(fch::Laplace<3>* laplace, const Medium& medium);

    void test_pic_vol2(fch::Laplace<3>* laplace, const Medium& medium, const TetgenMesh& mesh);

    void test_pic_vol3(const TetgenMesh& mesh) const;

    /** Interpolate electric field and potential on a Medium atoms */
    void interpolate(const Medium &medium);

    /** Interpolate electric field and potential on a set of points */
    void interpolate(const int n_points, const double* x, const double* y, const double* z);

    /** Calculate the electric field for the stationary current and temperature solver */
    void transfer_elfield(fch::CurrentsAndHeatingStationary<3>* ch_solver);

    /** Calculate the electric field for the transient current and temperature solver */
    void transfer_elfield(fch::CurrentsAndHeating<3>& ch_solver);

    /** Return electric field in i-th interpolation point */
    Vec3 get_elfield(const int i) const;

    /** Return electric field norm in i-th interpolation point */
    double get_elfield_norm(const int i) const;

    /** Return electric potential in i-th interpolation point */
    double get_potential(const int i) const;

    /** Compare the analytical and calculated field enhancement.
     * The check is disabled if lower and upper limits are the same. */
    bool check_limits(const vector<Solution>* solutions=NULL) const;

    /** Set parameters to calculate analytical solution */
    void set_check_params(const double E0, const double limit_min, const double limit_max,
            const double radius1, const double radius2=-1);

    /** Export calculated electic field distribution to HOLMOD */
    int export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    /** Interpolate electric field on set of points using the solution on surface mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    int interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric field on set of points using the solution on volumetric mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    int interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Interpolate electric potential on set of points using the solution on tetrahedral mesh nodes
     * @return  index of first point outside the mesh; index == -1 means all the points were inside the mesh */
    int interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
            double* phi, int* flag);

private:
    /** Data needed for comparing numerical solution with analytical one */
    double E0;                      ///< Long-range electric field strength
    double radius1;                 ///< Minor semi-axis of ellipse
    double radius2;                 ///< Major semi-axis of ellipse

    /** Return analytical electric field value for i-th point near the hemisphere */
    Vec3 get_analyt_field(const int i, const Point3& origin) const;

    /** Return analytical potential value for i-th point near the hemisphere */
    double get_analyt_potential(const int i, const Point3& origin) const;

    /** Get analytical field enhancement for hemi-ellipsoid on infinite surface */
    double get_analyt_enhancement() const;

    /** Interpolate electric field for heating module */
    void interpolate(vector<double>& elfields, const vector<dealii::Point<3>>& nodes);
};

/** Class to interpolate current densities and temperatures */
class HeatReader: public SolutionReader {
public:

    HeatReader(Interpolator* i);

    /** Interpolate solution on medium atoms */
    void interpolate(const Medium &medium);

    /** Linearly interpolate currents and temperatures in the bulk.
     *  In case of empty interpolator, constant values are stored. */
    void interpolate(fch::CurrentsAndHeating<3>& ch_solver);

    /** Export interpolated temperature */
    int export_temperature(const int n_atoms, double* T);

    Vec3 get_rho(const int i) const;

    double get_rho_norm(const int i) const;

    double get_temperature(const int i) const;

private:
};

/** Class to calculate field emission effects with GETELEC */
class EmissionReader: public SolutionReader {
public:

    EmissionReader(const FieldReader& fields, const HeatReader& heat, Interpolator* i);

    /** Calculates the emission currents and Nottingham heat distributions, including a rough
     * estimation of the space charge effects.
     * @param ch_solver heat solver object where J and Nottingham BCs will be written
     * @param conf Emission configuration parameters struct
     * @param Vappl Applied voltage (required for space charge calculations)
     */
    void calc_emission(const Config::Emission &conf, double Vappl = -1);

    void export_emission(fch::CurrentsAndHeating<3>& ch_solver);

    /** Injects electron SPs at the surface faces, depending on the current and the timestep */
    void inject_electrons(double delta_t, double Wsp, vector<Point3> &pos,
            vector<Vec3> &efield, vector<int> &cells);

    /** Initialises class data */
    void initialize(const TetgenMesh* m);

    double get_multiplier() const { return global_data.multiplier; }

    void set_multiplier(double _multiplier) { global_data.multiplier = _multiplier; }

private:
    /** Prepares the line inputed to GETELEC.
     *
     * @param point      Starting point of the line
     * @param direction  Direction of the line
     * @param rmax       Maximum distance that the line extends
     */
    void emission_line(const Point3& point, const Vec3& direction, const double rmax);

    /** Calculates representative quantities Jrep and Frep for space charge calculations
     * (See https://arxiv.org/abs/1710.00050 for definitions)
     */
    void calc_representative();

    /**
     * Calculates electron emission distribution for a given configuration (
     * @param workfunction  Input work function.
     */
    void emission_cycle(double workfunction, bool blunt  = false);

    static constexpr double electrons_per_fs = 6.2415e3; ///< definition of 1 ampere
    static constexpr double angstrom_per_nm = 10.0;
    static constexpr double nm2_per_angstrom2 = 0.01;
    static constexpr int n_lines = 32; ///< Number of points in the line for GETELEC

    const FieldReader& fields;    ///< Object containing the field on centroids of hex interface faces.
    const HeatReader& heat;       ///< Object containing the temperature on centroids of hexahedral faces.
    const TetgenMesh* mesh;     ///< Object containing information on the mesh.

    vector<double> current_densities;    ///< Vector containing the emitted current density on the interface faces [in Amps/A^2].
    vector<double> nottingham; ///< Same as current_densities for nottingham heat deposition [in W/A^2]
    vector<double> currents;    ///< Current flux for every face (current_densities * face_areas) [in Amps]

    vector<double> rline;   ///< Line distance from the face centroid (passed into GETELEC)
    vector<double> Vline;   ///< Potential on the straight line (complements rline)

    struct EmGlobalData {
        double multiplier;      ///< Multiplier for the field for Space Charge.
        double Jmax;    ///< Maximum current density of the emitter [in amps/A^2]
        double Fmax = 0.;    ///< Maximum local field on the emitter [V/A]
        double Frep = 0.;    ///< Representative local field (used for space charge equation) [V/A]
        double Jrep = 0.;    ///< Representative current deinsity for space charge. [amps/A^2]
        double I_tot = 0;   ///< Total current running through the surface [in Amps]
        double I_fwhm = 0;
    } global_data;
};

/** Class to calculate charges from electric field */
class ChargeReader: public SolutionReader {
public:

    ChargeReader(Interpolator* i);

    /** Calculate charge on the triangular faces using direct solution on the face centroid */
    void calc_charges(const TetgenMesh& mesh, const double E0);

    /** Remove the atoms and their solutions outside the MD box */
    void clean(const Medium::Sizes& sizes, const double latconst);

    /** Set parameters to calculate analytical solution */
    void set_check_params(const double Q_tot, const double limit_min, const double limit_max);

    /** Check whether charge is conserved within specified limits.
     * The check is disabled if lower and upper limits are the same. */
    bool check_limits(const vector<Solution>* solutions=NULL) const;

    /** Get electric field on the centroid of i-th triangle */
    Vec3 get_elfield(const int i) const;

    /** Get area of i-th triangle */
    double get_area(const int i) const;

    /** Get charge of i-th triangle */
    double get_charge(const int i) const;

private:
    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]
    double Q_tot;  ///< Analytical total charge
};

/** Class to calculate forces from charges and electric fields */
class ForceReader: public SolutionReader {
public:

    ForceReader(Interpolator* i);

    /** Calculate forces from atomic electric fields and face charges */
    void distribute_charges(const FieldReader &fields, const ChargeReader& faces,
        const double r_cut, const double smooth_factor);

    void calc_lorentz(const FieldReader &fields);

    /** Build Voronoi cells around the atoms in the region of interest */
    int calc_voronois(VoronoiMesh& mesh, const vector<int>& atom2face,
            const double radius, const double latconst, const string& mesh_quality);

    /** Calculate atomistic charges and Lorentz forces by using the Voronoi cells */
    void calc_charge_and_lorentz(const VoronoiMesh& mesh, const FieldReader& fields);

    /** Using the previously found surface charge, calculate Coulomb forces between atoms */
    void calc_coulomb(const double r_cut);

    /** Export the induced charge and force on imported atoms
     * @param n_atoms  number of first atoms the data will be exported
     * @param xq       charge and force in PARCAS format (xq[0] = q1, xq[1] = Fx1, xq[2] = Fy1, xq[3] = Fz1, xq[4] = q2, xq[5] = Fx2 etc)
     */
    int export_charge_and_force(const int n_atoms, double* xq) const;

    /** Export Laplace + Coulomb force and pair potential on imported atoms
     * @param n_atoms  number of first atoms the data will be exported
     * @param xnp      forces in PARCAS format & units (xnp[0] = Fx1, xnp[1] = Fy1, xnp[2] = Fz1, xnp[3] = Fx2 etc)
     * @param Epair    potential energy per atom
     * @param Vpair    total potential energy of atoms. Pot. due to Coloumb forces are added here. NOTE: Lorentz is missing!
     */
    int export_force_and_pairpot(const int n_atoms, double* xnp, double* Epair, double* Vpair) const;

    /** Return the force that is applied to i-th atom */
    Vec3 get_force(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].vector;
    }

    /** Return pair potential of i-th atom */
    double get_pairpot(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].norm;
    }

    /** Return the charge of i-th atom */
    double get_charge(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].scalar;
    }

private:
    static constexpr double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]
    static constexpr double force_factor = 0.5;  ///< force_factor = force / (charge * elfield)
    static constexpr double couloumb_constant = 14.399645; ///< force factor in Couloumb's law [V*A/e], == 1 / (4*pi*eps0)

    /** Screening factor for Coulomb force.
     * For details see Djurabekova et al, 2011, Physical Review E, 83(2), p.026704 */
    static constexpr double q_screen = 0.6809;

    /** Remove cells with too big faces*/
    void clean_voro_faces(VoronoiMesh& mesh);

    /** Separate cylindrical region from substrate region */
    int get_nanotip(Medium& nanotip, const double radius);
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
