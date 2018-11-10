/*
 * EmissionReader.h
 *
 *  Created on: 27.9.2018
 *      Author: kyritsak
 */
#ifndef EMISSIONREADER_H_
#define EMISSIONREADER_H_

#include "SolutionReader.h"
#include "getelec.h"


using namespace std;
namespace femocs {

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
    void calc_emission(const Config::Emission &conf, double Veff_SC = -1,
            bool update_eff_region = false);

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

    vector<double>* get_current_densities() {
        return &current_densities;
    }

    vector<double>* get_nottingham() {
        return &nottingham;
    }

    /**
     * Get the statistic data for a dat file
     * @param first_line
     * @return
     */
    string get_stats(const bool first_line) const;

    /** Compose entry to dat file */
    string get_global_data(const bool first_line) const;

    struct EmGlobalData {
        double multiplier=1.;   ///< Multiplier for the field for Space Charge.
        double theta=1.;       ///< correction multiplier for Space Charge
        double sfactor=1.;     ///< factor that rescales the field (accounts for changed applied V or F)
        double Jmax;         ///< Maximum current density of the emitter [in amps/A^2]
        double Fmax = 0.;    ///< Maximum local field on the emitter [V/A]
        double Frep = 0.;    ///< Representative local field (used for space charge equation) [V/A]
        double Jrep = 0.;    ///< Representative current deinsity for space charge. [amps/A^2]
        double I_tot = 0;    ///< Total current running through the surface [in Amps]
        double I_eff = 0;   ///< Total current within the effective area
        double area;        ///< total area of emitting region
    } global_data;


    struct EmStats{
        int N_calls;        ///< Counter keeping the last N_calls

        vector<double> I_tot; ///< List of all the I_tot for the last N_calls (useful for convergence check)
        double Itot_mean = 0;  ///< Mean current of the last Ilist;
        double Itot_std = 0;  ///< STD of the last Ilist;

        vector<double> Jrep; ///< List of all the Jrep for the last N_calls (useful for convergence check)
        double Jrep_mean = 0;  ///< Mean of the last Jrep list;
        double Jrep_std = 0;  ///< STD of the last Jrep list;

        vector<double> Frep; ///< List of all the Frep for the last N_calls (useful for convergence check)
        double Frep_mean = 0;  ///< Mean of the last Frep list;
        double Frep_std = 0;  ///< STD of the last Frep list;

        vector<double> Jmax; ///< List of all the Jmax for the last N_calls (useful for convergence check)
        double Jmax_mean = 0;  ///< Mean of the last Jmax list;
        double Jmax_std = 0;  ///< STD of the last Jmax list;

        vector<double> Fmax; ///< List of all the Fmax for the last N_calls (useful for convergence check)
        double Fmax_mean = 0;  ///< Mean of the last Fmax list;
        double Fmax_std = 0;  ///< STD of the last Fmax list;

    } stats;

private:
    /** Prepares the line inputed to GETELEC.
     *
     * @param point      Starting point of the line
     * @param direction  Direction of the line
     * @param rmax       Maximum distance that the line extends
     */
    void emission_line(const Point3& point, const Vec3& direction, const double rmax);

    /** Calculates all the global values */
    void calculate_globals();

    /** Compose console output of occured errors */
    string get_error_codes(vector<int> &errors) const;

    /**
     * Calculates the effective emission area (assigns flag to each surface face
     * @param threshold minimum value (fraction of maximum) to be considered effective area
     * @param mode "field" or "current", whether the criterion is on the field or the current
     */
    void calc_effective_region(double threshold, string mode);

    /** Compose entry to xyz or movie file */
    string get_data_string(const int i) const;

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
    vector<bool> is_effective;  ///< effective emission area
    vector<double> rline;   ///< Line distance from the face centroid (passed into GETELEC)
    vector<double> Vline;   ///< Potential on the straight line (complements rline)
    vector<double> thetas_SC; ///< local field reduction factor due to SC

    friend class Pic<3>;   // for convenience, allow Pic-class to access private data
};

}

#endif /* EMISSIONREADER_H_ */
