/*
 * physical_quantities.h
 *
 *  Created on: Apr 30, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_PHYSICAL_QUANTITIES_H_
#define INCLUDE_PHYSICAL_QUANTITIES_H_

#include <vector>
#include <string>
#include <utility>
#include "Config.h"

namespace fch {

/** @brief Evaluation tools for emission and conductivities based on tabulated data.
 * Resisitivity data is from literature, emission and nottingham data is pre-calculated
 * and tabulated based on WKB approximation, Schottky-Nordheim barrier,
 * and free electron model.
 * See Kristjan Eimre MSc thesis for details.
 */
class PhysicalQuantities {
public:
    /**
     * The default constructor initializes the physical values to hardcoded ones.
     * Can be overwritten by loading the data from a file.
     */
    PhysicalQuantities(const femocs::Config::Heating& config);

    /**
     * Load emission current data from file
     * @return true if successful, false otherwise
     */
    bool load_emission_data(std::string filepath);

    /**
     * Load nottingham delta energy data from file
     * @return true if successful, false otherwise
     */
    bool load_nottingham_data(std::string filepath);

    /**
     * Load copper resistivity data from file
     * @return true if successful, false otherwise
     */
    bool load_resistivity_data(std::string filepath);

    /**
     * Evaluates the emission current density J
     * @param field electric field in (GV/m)
     * @param temperature Temperature in (K)
     * @return Emission current density in (A/ang^2)
     */
    double emission_current(double field, double temperature);

    /**
     * Evaluates the nottingham delta energy <DE> (eV)
     * @param field electric field in (GV/m)
     * @param temperature Temperature in (K)
     * @return nottingham delta energy in (eV)
     */
    double nottingham_de(double field, double temperature);

    /**
     * Evaluates the electrical resistivity rho
     * @param temperature T in (K)
     * @return resistivity in (Ohm*ang)
     */
    double evaluate_resistivity(double temperature) const;

    double evaluate_resistivity_derivative(double temperature);

    /**
     * electrical conductivity sigma in (1/(Ohm*ang))
     */
    double sigma(double temperature) const;

    /**
     * electrical conductivity derivative
     */
    double dsigma(double temperature);

    /**
     * thermal conductivity in (W/(ang*K))
     */
    double kappa(double temperature);

    /**
     * thermal conductivity derivative
     */
    double dkappa(double temperature);

    /**
     * Outputs sigma, kappa, res (and d-s) and emission currents to files in specified path
     * NB: Slow!!!
     */
    void output_to_files();

    /**
     * Method to read the tabulated resistivity data from file.
     * If file not found, hard coded data are used.
     */
    void initialize_with_hc_data();

private:
    const femocs::Config::Heating& config;
    /**
     * A data structure to hold the uniform grid information,
     * used in bilinear interpolation
     * Holds 2d information in a 1d vector: to access element v[i][j], use v[i*ynum+j]
     */
    struct InterpolationGrid {
        std::vector<double> v;
        double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
        unsigned xnum = 0, ynum = 0;
    };
    InterpolationGrid emission_grid;
    InterpolationGrid nottingham_grid;

    std::vector<std::pair<double, double> > resistivity_data;

    /**
     * 1d linear interpolation with constant extrapolation using binary search
     */
    double linear_interp(double x,
            std::vector<std::pair<double, double>> data) const;

    /**
     * 1d linear interpolation of the derivative with constant extrapolation using binary search
     * derivative is approximated with central differences (one sided at ends)
     */
    double deriv_linear_interp(double x,
            std::vector<std::pair<double, double>> data);

    double evaluate_derivative(std::vector<std::pair<double, double>> &data,
            std::vector<std::pair<double, double>>::iterator it);

    /**
     * 2d bilinear interpolation with constant extrapolation
     * NB: Assumes uniform grid
     */
    double bilinear_interp(double x, double y,
            const InterpolationGrid &grid_data);

    bool load_spreadsheet_grid_data(std::string filepath,
            InterpolationGrid &grid);

    bool load_compact_grid_data(std::string filepath, InterpolationGrid &grid);

    /**
     * Hardcoded data
     */
    static const std::vector<std::pair<double, double>> hc_resistivity_data;
    static const std::vector<double> hc_emission_current_data;
    static const std::vector<double> hc_nottingham_data;

};

} // namespace fch

#endif /* INCLUDE_PHYSICAL_QUANTITIES_H_ */
