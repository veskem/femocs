/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include "Primitives.h"
#include "Surface.h"
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
    SolutionReader(Interpolator* i, const string& vec_lab, const string& scalar1_lab, const string& scalar2_lab);

    /** Interpolate solution using already available data about which atom is connected with which cell */
    virtual void calc_interpolation();

    /** Map atoms to cells and interpolate solution on the system atoms */
    void calc_full_interpolation();

    /** Reserve memory for data */
    void reserve(const int n_nodes);

    /** Print statistics about interpolated solution */
    void print_statistics();

    /** Get pointer to interpolation vector */
    vector<Solution>* get_interpolations() { return &interpolation; }

    /** Get i-th Solution */
    Solution get_interpolation(const int i) const {
        require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + d2s(i));
        return interpolation[i];
    }

    /** Set i-th Solution */
    void set_interpolation(const int i, const Solution& s) {
        require(i >= 0 && i < static_cast<int>(interpolation.size()), "Index out of bounds: " + d2s(i));
        interpolation[i] = s;
    }

    /** Set interpolation preferences */
    void set_preferences(const bool _srt, const int _dim, const int _rank, const bool _centroid=false) {
        require((_dim == 2 || _dim == 3), "Invalid interpolation dimension: " + d2s(_dim));
        require((_rank == 1 || _rank == 2 || _rank == 3), "Invalid interpolation rank: " + d2s(_rank));
        sort_atoms = _srt;
        interp_centroids = _centroid;
        dim = _dim;
        rank = _rank;
    }

    /** Alter the pointer to interpolator */
    void set_interpolator(Interpolator* i) { interpolator = i; }

    /** Interpolate solution on the surface mesh centroids of the FEM solver */
    void interpolate(const DealSolver<3>& solver);

    /** Interpolate solution on a set of given points */
    void interpolate(const int n_points, const double* x, const double* y, const double* z);

    /** Determine whether given data is included in SolutionReader */
    int contains(const string& data_label) const;

    /** Determine whether ref_label or its lower case copy matches with label1 or label2 */
    int contains(const string& ref_label, const string& label1, const string& label2="") const;

    /** General function to export desired component of calculated data */
    int export_results(const int n_points, const string &data_type, double* data) const;

    /** General function to first perform interpolation
     * and then export desired component of calculated data */
    int interpolate_results(const int n_points, const string &data_type, const double* x,
            const double* y, const double* z, double* data);

protected:
    const string vec_label;       ///< label for vector data
    const string scalar1_label;   ///< label for first scalar data
    const string scalar2_label;   ///< label for second scalar data
    double limit_min;             ///< minimum value of accepted comparison value
    double limit_max;             ///< maximum value of accepted comparison value
    bool sort_atoms;              ///< sort atoms along Hilbert curve to make interpolation faster
    bool interp_centroids;        ///< points are centroids of cells
    int dim;                      ///< location of interpolation; 2-surface, 3-space
    int rank;                     ///< interpolation rank; 1-linear, 2-quadratic
    bool atoms_mapped_to_cells;   ///< flag indicating whether fast interpolation can be performed or not

    Interpolator* interpolator;    ///< pointer to interpolator
    vector<Solution> interpolation;       ///< interpolated data

    /** Output atom data in .xyz format */
    void write_xyz(ofstream &outfile) const;

    /** Output information about data vectors into .vtk file. */
    void write_vtk_point_data(ofstream& out) const;

    /** Find the cell that surrounds i-th atom and interpolate the solution for it.
     * @param cell   initial guess for the cell that might surround the point
     * @return cell index that was found to surround the point */
    int locate_interpolate(const int i, int cell);

    /** Find the centroid that is closest to the i-th point and interpolate the solution for it */
    int locate_interp_centroid(const int i, int cell);

    /** Interpolate the solution for i-th point */
    void interp_solution(const int i);

    /** Interpolate the solution for i-th centroidal point */
    void interp_centroid(const int i);

    /** Sort atoms and interpolation by the atom ID */
    void restore_sorting();

    /** Export vector component of solution */
    void export_vec(const int n_points, double* data, bool append) const;

    /** Export norm component of solution */
    void export_norm(const int n_points, double* data, bool append) const;

    /** Export scalar component of solution */
    void export_scalar(const int n_points, double* data, bool append) const;

    /** Clear data vector */
    void clear_data(const int n_data, double* data) const;
};

/** Class to extract solution from DealII calculations */
class FieldReader: public SolutionReader {
public:

    FieldReader(Interpolator* i);

    /** Interpolate solution on all Medium atoms */
    using SolutionReader::interpolate;
    void interpolate(const Surface &surface);

    /** Compare the analytical and calculated field enhancement.
     * The check is disabled if lower and upper limits are the same. */
    bool check_limits(const vector<Solution>* solutions=NULL, bool verbose=true);

    /** Set parameters to calculate analytical solution */
    void set_check_params(const Config& conf, double radius, double tip_height, double box_height);

    /** In addition to regular interpolation, pre-calculate also field norms */
    void calc_interpolation();

    /** Return electric field in i-th interpolation point */
    Vec3 get_elfield(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].vector;
    }

    /** Return electric field norm in i-th interpolation point */
    double get_elfield_norm(const int i) const {
        require(i >= 0 && i < field_norm.size(), "Invalid index: " + d2s(i));
        return field_norm[i];
    }

    /** Return charge density in i-th interpolation point */
    double get_charge_dens(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar1;
    }

    /** Return electric potential in i-th interpolation point */
    double get_potential(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar2;
    }

    /** Find the maximum field norm either from the provided of internal solution vector */
    double max_field(const vector<Solution>* solutions=NULL) const;

    /** Return measured-"analytical" field enhancement ratio that was calculated in check_limits */
    double get_beta() const { return beta; }

private:
    vector<double> field_norm;  ///< pre-computed electric field norms

    /** Data needed for comparing numerical solution with analytical one */
    double E0;         ///< Long-range electric field strength
    double radius1;    ///< Minor semi-axis of ellipse
    double radius2;    ///< Major semi-axis of ellipse
    double beta;       ///< Measured / "analytical" field enhancement
    double E_max;      ///< Maximum electric field norm value

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

    /** Interpolate solution on non-fixed AtomReader atoms */
    using SolutionReader::interpolate;
    void interpolate(const AtomReader &reader);

    /** Interpolate and export solution on the mesh DOFs of the FEM solver */
    void interpolate_dofs(CurrentHeatSolver<3>& solver, const TetgenMesh* mesh);

    /** Compute data that Berendsen thermostat requires for re-using old solution */
    void precalc_berendsen(bool update_locations);

    /** Apply Berendsen thermostat for atomistic velocities */
    int scale_berendsen(double* x1, const int n_atoms, const vector<Vec3>& velocities, const Config& conf);

    /** Export kinetic energy that was added due to Berendsen velocity scaling */
    int export_kin_energy(const string &data_type, double* energy) const;

    /** Return current density in i-th interpolation point */
    Vec3 get_rho(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].vector;
    }

    /** Return magnitude of current density in i-th interpolation point */
    double get_potential(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar1;
    }

    /** Return temperature in i-th interpolation point */
    double get_temperature(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar2;
    }

private:
    static constexpr double kB = 8.6173303e-5; ///< Boltzmann constant [eV/K]
    static constexpr double amu = 103.642696;  ///< 1 amu in eV*(fs/A)^2
    static constexpr double heat_factor = amu / (2*1.5*kB);  ///< Factor to transfer 2*kinetic energy to temperature

    vector<vector<int>> tet2atoms;
    vector<double> fem_temp;
    vector<double> temperatures;
    vector<double> lambdas;
    double kin_energy;              ///< Kinetic energy added due to Berendsen scaling [eV]

    /** Output atom data in .xyz format */
    void write_xyz(ofstream &outfile) const;

    /** Function to increase spatial ordering of mesh nodes */
    void sort_spatial();
    void sort_spatial(const TetgenMesh* mesh);

    /** Calculate temperature scaling factors for Berendsen thermostat */
    void calc_lambdas(const vector<Vec3>& velocities, const Config& conf);
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
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].vector;
    }

    /** Get area of i-th triangle */
    double get_area(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar1;
    }

    /** Get charge of i-th triangle */
    double get_charge(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar2;
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
    int calc_voronois(VoronoiMesh& mesh, const Coarseners& coarseners, const Config::Geometry& conf, const string& mesh_quality);

    /** Calculate Lorentz forces from known electric fields and charges */
    void recalc_lorentz(const FieldReader &fields);

    /** Calculate atomistic charges and Lorentz forces by using the Voronoi cells */
    void calc_charge_and_lorentz(const VoronoiMesh& mesh, const FieldReader& fields);

    /** Using the previously found surface charge, calculate Coulomb forces between atoms */
    void calc_coulomb(const double r_cut);

    /** Calculate Coulomb forces from known charges and atom neighbouring */
    void recalc_coulomb();

    /** Export Parcas-specific data */
    int export_parcas(const int n_points, const string &data_type, const Vec3& si2parcas, double* data) const;

    /** Export sum of pair-potential of Coulomb interaction */
    int export_pot_energy(const int n_points, const string &data_type, double* energy) const;

    /** Return the force that is applied to i-th atom */
    Vec3 get_force(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].vector;
    }

    /** Return pair potential of i-th atom */
    double get_pairpot(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar1;
    }

    /** Return the charge of i-th atom */
    double get_charge(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return interpolation[i].scalar2;
    }

private:
    static constexpr double eps0 = 0.0055263494; ///< vacuum permittivity [e/V*A]
    static constexpr double force_factor = 0.5;  ///< force_factor = force / (charge * elfield)
    static constexpr double couloumb_constant = 14.399645; ///< force factor in Couloumb's law [V*A/e], == 1 / (4*pi*eps0)

    /** Screening factor for Coulomb force.
     * For details see Djurabekova et al, 2011, Physical Review E, 83(2), p.026704 */
    static constexpr double q_screen = 0.6809;

    vector<int> nborlist;   ///< diagonal neighbour list to speed up Coulomb recalculation

    /** Remove cells with too big faces*/
    void clean_voro_faces(VoronoiMesh& mesh);

    /** Separate region-of-interest from substrate region */
    int get_nanotip(Medium& nanotip, const Coarseners& coarseners);
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
