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

class PhysicalQuantities {
public:

	bool load_emission_data(std::string filepath);
	bool load_resistivity_data(std::string filepath);

	/**
	 * Evaluates the emission current
	 * @param field electric field in (GV/m)
	 * @param temperature Temperature in (K)
	 * @return Emission current density in (A/nm^2)
	 */
	double evaluate_current(double field, double temperature);

	/**
	 * Evaluates the electrical conductivity
	 * @param temperature Temperature in (K)
	 * @return resistivity in (Ohm*m)
	 */
	double evaluate_resistivity(double temperature);

	double evaluate_resistivity_derivative(double temperature);

	/**
	 * electrical conductivity in (1/(Ohm*nm))
	 */
	double sigma(double temperature);

	/**
	 * electrical conductivity derivative
	 */
	double dsigma(double temperature);

	/**
	 * thermal conductivity in (W/(nm*K))
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

private:
	/**
	 * A data structure to hold the uniform grid information,
	 * used in bilinear interpolation
	 */
	struct InterpolationGrid {
		std::vector<std::vector<double> > grid;
		double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
	};
	InterpolationGrid emission_grid;

	std::vector<std::pair<double, double> > resistivity_data;

	/**
	 * 1d linear interpolation with constant extrapolation using binary search
	 */
	double linear_interp(double x, std::vector<std::pair<double, double>> data);

	/**
	 * 1d linear interpolation of the derivative with constant extrapolation using binary search
	 * derivative is approximated with central differences (one sided at ends)
	 */
	double deriv_linear_interp(double x, std::vector<std::pair<double, double>> data);

	double evaluate_derivative(std::vector<std::pair<double, double>> &data,
			std::vector<std::pair<double, double>>::iterator it);

	/**
	 * 2d bilinear interpolation with constant extrapolation
	 * NB: Assumes uniform grid
	 */
	double bilinear_interp(double x, double y, InterpolationGrid grid_data);

};

#endif /* INCLUDE_PHYSICAL_QUANTITIES_H_ */
