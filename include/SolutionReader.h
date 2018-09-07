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
#include "AtomReader.h"
#include "TetgenCells.h"
#include "TetgenMesh.h"
#include "VoronoiMesh.h"
#include "Interpolator.h"
#include "Config.h"
#include "CurrentHeatSolver.h"
#include "DealSolver.h"

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

    /** Get pointer to interpolation vector */
    vector<Solution>* get_interpolations() { return &interpolation; }

    /** Get i-th Solution */
    Solution get_interpolation(const int i) const {
        require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + to_string(i));
        return interpolation[i];
    }

    /** Set i-th Solution */
    void set_interpolation(const int i, const Solution& s) {
        require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + to_string(i));
        interpolation[i] = s;
    }

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

    /** Store the surface mesh centroids of the FEM solver */
    void store_points(const DealSolver<3>& solver);

    /** Interpolate solution on the surface mesh centroids of the FEM solver */
    void interpolate(const DealSolver<3>& solver);

    /** Interpolate solution on a set of given points */
    void interpolate(const int n_points, const double* x, const double* y, const double* z);

    /** Interpolate solution on all Medium atoms */
    void interpolate(const Medium &medium);

    /** Interpolate solution on non-fixed AtomReader atoms */
    void interpolate(const AtomReader &reader);

    /** Determine whether given data is included in SolutionReader */
    int contains(const string& data_label) const;

    /** General function to export desired component of calculated data */
    int export_results(const int n_points, const string &data_type, double* data) const;

    /** General function to first perform interpolation
     * and then export desired component of calculated data */
    int interpolate_results(const int n_points, const string &data_type, const double* x,
            const double* y, const double* z, double* data);

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

    /** Find the cell that surrounds i-th atom and interpolate the solution for it.
     * @param cell   initial guess for the cell that might surround the point
     * @return cell index that was found to surround the point */
    int update_interpolation(const int i, int cell);
};

/** Class to extract solution from DealII calculations */
class FieldReader: public SolutionReader {
public:

    FieldReader(Interpolator* i);

    /** Compare the analytical and calculated field enhancement.
     * The check is disabled if lower and upper limits are the same. */
    bool check_limits(const vector<Solution>* solutions=NULL) const;

    /** Find the maximum field norm from the solution vector */
    double calc_max_field(const vector<Solution>* solutions=NULL) const;

    /** Set parameters to calculate analytical solution */
    void set_check_params(const Config& conf, double radius, double tip_height, double box_height);

    /** Return electric field in i-th interpolation point */
    Vec3 get_elfield(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].vector;
    }

    /** Return electric field norm in i-th interpolation point */
    double get_elfield_norm(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].norm;
    }

    /** Return electric potential in i-th interpolation point */
    double get_potential(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].scalar;
    }

private:
    /** Data needed for comparing numerical solution with analytical one */
    double E0;         ///< Long-range electric field strength
    double radius1;    ///< Minor semi-axis of ellipse
    double radius2;    ///< Major semi-axis of ellipse

    /** Return analytical electric field value for i-th point near the hemisphere */
    Vec3 get_analyt_field(const int i, const Point3& origin) const;

    /** Return analytical potential value for i-th point near the hemisphere */
    double get_analyt_potential(const int i, const Point3& origin) const;

    /** Get analytical field enhancement for hemi-ellipsoid on infinite surface */
    double get_analyt_enhancement() const;
};

/** Class to interpolate current densities and temperatures */
class HeatReader: public SolutionReader {
public:

    HeatReader(Interpolator* i);

    void locate_atoms(const Medium &medium);

    /** Compute data that Berendsen thermostat requires for re-using old solution */
    void precalc_berendsen_long();

    /** Apply Berendsen thermostat for individual atoms */
    int scale_berendsen_short(double* x1, const int n_atoms, const Vec3& parcas2si);

    /** Apply Berendsen thermostat for atoms within a tetrahedron */
    int scale_berendsen_long(double* x1, const int n_atoms, const Vec3& parcas2si);

    /** Store velocity scaling constants */
    void set_params(const Config& conf) {
        data.tau = conf.heating.tau;
        data.md_timestep = conf.behaviour.timestep_fs;
        data.time_unit = 10.1805*sqrt(conf.behaviour.mass);
    }

    /** Return current density in i-th interpolation point */
    Vec3 get_rho(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].vector;
    }

    /** Return magnitude of current density in i-th interpolation point */
    double get_rho_norm(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].norm;
    }

    /** Return temperature in i-th interpolation point */
    double get_temperature(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].scalar;
    }

private:
    static constexpr double kB = 8.6173324e-5; ///< Boltzmann constant [eV/K]
    static constexpr double heat_factor = 1.0 / (2*1.5*kB);  ///< Factor to transfer 2*kinetic energy to temperature

    vector<vector<int>> tet2atoms;
    vector<double> fem_temp;

    struct Data {
        double tau;          ///< Time constant in Berendsen scaling [fs]
        double md_timestep;  ///< MD time step [fs]
        double time_unit;    ///< the conversion factor of Parcas internal units to fs.
    } data;

    /** Transfer velocities from Parcas units to fm / fs */
    void calc_SI_velocities(vector<Vec3>& velocities, const int n_atoms, const Vec3& parcas2si, double* x1);

    /** Calculate scaling factor for Berendsen thermostat */
    double calc_lambda(const double T_start, const double T_end) const;
};

// forward declaration of Pic for declaring it as a friend
template<int dim> class Pic;

/** Class to calculate field emission effects with GETELEC */
class EmissionReader: public SolutionReader {
public:

    EmissionReader(const FieldReader *fields, const HeatReader *heat,
            const PoissonSolver<3> *poisson, Interpolator* i);

    /** Calculates the emission currents and Nottingham heat distributions, including a rough
     * estimation of the space charge effects.
     * @param ch_solver heat solver object where J and Nottingham BCs will be written
     * @param conf Emission configuration parameters struct
     * @param Veff_SC (effective applied voltage for space charge calculations)
     */
    void calc_emission(const Config::Emission &conf, double Veff_SC = -1);

    void export_emission(CurrentHeatSolver<3>& ch_solver);

    /**
     * Calculates the mean and the standard deviation of the total current for the last N_calls
     * @param std the returned standard deviation
     * @return
     */
    void calc_global_stats(void);

    /** Initialises class data */
    void initialize(const TetgenMesh* m, bool reinit = true);

    void set_sfactor(double factor){
        global_data.multiplier *= factor / global_data.sfactor;
        global_data.sfactor = factor;
    }

    struct EmGlobalData {
        double multiplier=1.;   ///< Multiplier for the field for Space Charge.
        double theta=1.;       ///< correction multiplier for Space Charge
        double sfactor=1.;     ///< factor that rescales the field (accounts for changed applied V or F)
        double Jmax;         ///< Maximum current density of the emitter [in amps/A^2]
        double Fmax = 0.;    ///< Maximum local field on the emitter [V/A]
        double Frep = 0.;    ///< Representative local field (used for space charge equation) [V/A]
        double Jrep = 0.;    ///< Representative current deinsity for space charge. [amps/A^2]
        double I_tot = 0;    ///< Total current running through the surface [in Amps]
        double I_fwhm = 0;   ///< Total current within the FWHM area
        int N_calls;        ///< Counter keeping the last N_calls
        double area;        ///< total area of emitting region
        vector<double> Ilist; ///< List of all the I_tot for the last N_calls (useful for convergence check)
        double I_mean = 0;  ///< Mean current of the last Ilist;
        double I_std = 0;  ///< STD of the last Ilist;
    } global_data;

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
    void emission_cycle(double workfunction, bool blunt  = false,
            bool cold = false, double Vappl = 0);

    /** Compose entry to xyz or movie file */
    string get_data_string(const int i) const;

    /** Compose entry to dat file */
    string get_global_data(const bool first_line) const;

    static constexpr double angstrom_per_nm = 10.0;
    static constexpr double nm2_per_angstrom2 = 0.01;
    static constexpr int n_lines = 32; ///< Number of points in the line for GETELEC

    const FieldReader *fields;      ///< field on centroids of hex interface faces.
    const HeatReader *heat;         ///< temperature on centroids of hexahedral faces.
    const TetgenMesh *mesh;         ///< data & operations about the mesh.
    const PoissonSolver<3> *poisson; ///< Poisson solver information (field, potential etc)

    vector<double> current_densities;    ///< Vector containing the emitted current density on the interface faces [in Amps/A^2].
    vector<double> nottingham; ///< Same as current_densities for nottingham heat deposition [in W/A^2]
    vector<double> currents;    ///< Current flux for every face (current_densities * face_areas) [in Amps]
    vector<double> rline;   ///< Line distance from the face centroid (passed into GETELEC)
    vector<double> Vline;   ///< Potential on the straight line (complements rline)

    friend class Pic<3>;   // for convenience, allow Pic-class to access private data
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
    Vec3 get_elfield(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].vector;
    }

    /** Get area of i-th triangle */
    double get_area(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].norm;
    }

    /** Get charge of i-th triangle */
    double get_charge(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return interpolation[i].scalar;
    }

private:
    const double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]
    double Q_tot;  ///< Analytical total charge
};

/** Class to calculate forces from charges and electric fields */
class ForceReader: public SolutionReader {
public:

    ForceReader(Interpolator* i);

    /** Calculate forces from atomic electric fields and face charges */
    void distribute_charges(const FieldReader &fields, const ChargeReader& faces, const double smooth_factor);

    /** Build Voronoi cells around the atoms in the region of interest */
    int calc_voronois(VoronoiMesh& mesh, const vector<int>& atom2face,
            const double radius, const double latconst, const string& mesh_quality);

    /** Calculate Lorentz forces from atomistic electric fields and face charges */
    void calc_lorentz(const FieldReader &fields);

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

    int export_parcas(const int n_points, const string &data_type, const Vec3& si2parcas, double* data) const;

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


/**
 * Class for performing tests with shape functions
 */
class InterpolatorTester : public SolutionReader {
public:
    InterpolatorTester(Interpolator* i);

    void compare_shape_funs(PoissonSolver<3> &poisson, const Medium::Sizes &sizes);
    void compare_interpolators(PoissonSolver<3> &poisson, const Medium::Sizes &sizes);
    void test_corners(const TetgenMesh& mesh) const;
    void compare_space(const Medium::Sizes &sizes);
    void compare_surface(const Medium &medium);
    void perform_comparison(const string &file_name);

private:
};
} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
