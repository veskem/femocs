/*
 * physical_quantities.cc
 *
 *  Created on: Apr 30, 2016
 *      Author: kristjan
 */

#include "physical_quantities.h"

#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>  // fopen
#include <iostream>
#include <algorithm>

double PhysicalQuantities::evaluate_current(double field, double temperature) {
	return std::exp(bilinear_interp(std::log(field), temperature, emission_grid))*1.0e-18;
}

double PhysicalQuantities::evaluate_resistivity(double temperature) const {
	return linear_interp(temperature, resistivity_data);
}

double PhysicalQuantities::evaluate_resistivity_derivative(double temperature) {
	return deriv_linear_interp(temperature, resistivity_data);
}

double PhysicalQuantities::sigma(double temperature) const {
	if (temperature < 200) temperature = 200;
	if (temperature > 1400) temperature = 1400;
	double rho = evaluate_resistivity(temperature)*1.0e9;
	return 1.0 / rho;
}

double PhysicalQuantities::dsigma(double temperature) {
	if (temperature < 200) temperature = 200;
	if (temperature > 1400) temperature = 1400;
	double rho = evaluate_resistivity(temperature)*1.0e9;
	return -evaluate_resistivity_derivative(temperature)*1.0e9 / (rho * rho);
}

double PhysicalQuantities::kappa(double temperature) {
	if (temperature < 200) temperature = 200;
	if (temperature > 1400) temperature = 1400;
	double lorentz = 2.443e-8;
	return lorentz * temperature * sigma(temperature);
}

double PhysicalQuantities::dkappa(double temperature) {
	if (temperature < 200) temperature = 200;
	if (temperature > 1400) temperature = 1400;
	double lorentz = 2.443e-8;
	return lorentz * (sigma(temperature) + temperature * dsigma(temperature));
}

bool PhysicalQuantities::load_emission_data(std::string filepath) {
	std::ifstream infile(filepath);
	if (!infile) {
		std::cerr << "Couldn't open \"" << filepath << "\"\n";
		return false;
	}

	std::vector<double> row;
	std::string line;
	double last_x = 0, last_y = 0;
	bool first_line = 1;

	while (std::getline(infile, line)) {
		if (line[0] == '%')
			continue;
		std::istringstream stm(line);
		double x, y, z;
		stm >> x >> y >> z;

		if (first_line) {
			emission_grid.xmin = x;
			emission_grid.ymin = y;
			first_line = 0;
		} else if (last_x != x) {
			emission_grid.grid.push_back(row);
			row.clear();
		}
		row.push_back(z);
		last_x = x;
		last_y = y;
	}
	emission_grid.grid.push_back(row);
	emission_grid.xmax = last_x;
	emission_grid.ymax = last_y;

	infile.close();
	return true;
}

bool PhysicalQuantities::load_resistivity_data(std::string filepath) {
	std::ifstream infile(filepath);
	if (!infile) {
		std::cerr << "Couldn't open \"" << filepath << "\"\n";
		return false;
	}
	double x, y;
	while (infile >> x >> y) {
		resistivity_data.push_back(std::make_pair(x, y));
	}
	infile.close();
	return true;
}

double PhysicalQuantities::linear_interp(double x, std::vector<std::pair<double, double>> data) const {
	if (x <= data[0].first)
		return data[0].second;
	if (x >= data.back().first)
		return data.back().second;
	typedef std::pair<double, double> myPair;
	auto pair_comp = [](myPair lhs, myPair rhs) -> bool {return lhs.first < rhs.first;};
// use binary search for the position:
	auto it1 = std::lower_bound(data.begin(), data.end(), std::make_pair(x, 0.0), pair_comp);
	auto it2 = it1 - 1;
	return it2->second + (it1->second - it2->second) * (x - it2->first) / (it1->first - it2->first);
}

double PhysicalQuantities::evaluate_derivative(std::vector<std::pair<double, double>> &data,
		std::vector<std::pair<double, double>>::iterator it) {
	if (it == data.begin()) {
		return ((it + 1)->second - it->second) / ((it + 1)->first - it->first);
	} else if (it == data.end()-1) {
		return (it->second - (it - 1)->second) / (it->first - (it - 1)->first);
	}
	return ((it + 1)->second - (it - 1)->second) / ((it + 1)->first - (it - 1)->first);
}

void PhysicalQuantities::output_to_files() {

	// temperature for emission current evaluation
	double temperature = 500.0;
	std::string filepath = "./output/";

	std::cout << "Outputting physical quantities to \""+filepath+"\"" << std::endl;

	FILE *rho_file = fopen((filepath+"rho_file.txt").c_str(), "w");
	FILE *sigma_file = fopen((filepath+"sigma_file.txt").c_str(), "w");
	FILE *kappa_file = fopen((filepath+"kappa_file.txt").c_str(), "w");
	FILE *emission_file = fopen((filepath+"emission_file.txt").c_str(), "w");

	for (double t = 100.0; t < 1500.0; t+=5.0) {
		fprintf(rho_file, "%.5e %.5e %.5e\n", t,
				evaluate_resistivity(t), evaluate_resistivity_derivative(t));
		fprintf(sigma_file, "%.5e %.5e %.5e\n", t, sigma(t), dsigma(t));
		fprintf(kappa_file, "%.5e %.5e %.5e\n", t, kappa(t), dkappa(t));
	}

	for (double f = 0.01; f < 10.0; f+=0.01) {
		fprintf(emission_file, "%.5e %.5e %.16e\n", f, temperature, evaluate_current(f, temperature));
	}

	fclose(rho_file);
	fclose(sigma_file);
	fclose(kappa_file);
	fclose(emission_file);
}

/**
 * NB: Derivative is extrapolated by boundary values; out of bounds the real derivative should be zero!
 */
double PhysicalQuantities::deriv_linear_interp(double x, std::vector<std::pair<double, double>> data) {
	double eps = 1e-10;
	if (x <= data[0].first)
		x = data[0].first + eps;
	if (x >= data.back().first)
		x = data.back().first;
	typedef std::pair<double, double> myPair;
	auto pair_comp = [](myPair lhs, myPair rhs) -> bool {return lhs.first < rhs.first;};
// use binary search for the position:
	auto it = std::lower_bound(data.begin(), data.end(), std::make_pair(x, 0.0), pair_comp);
	auto itp = it - 1;

	double dit = evaluate_derivative(data, it);
	double ditp = evaluate_derivative(data, itp);

	return ditp + (dit - ditp) * (x - itp->first) / (it->first - itp->first);
}

// NB: This assumes uniform grid
double PhysicalQuantities::bilinear_interp(double x, double y, const InterpolationGrid &grid_data) {
//std::printf("%f, %f, %f\n", x, grid_data.xmin, grid_data.xmax);
//std::printf("%f, %f, %f\n", y, grid_data.ymin, grid_data.ymax);
	double eps = 1e-10;
	if (x <= grid_data.xmin)
		x = grid_data.xmin;
	if (x >= grid_data.xmax)
		x = grid_data.xmax - eps;
	if (y <= grid_data.ymin)
		y = grid_data.ymin;
	if (y >= grid_data.ymax)
		y = grid_data.ymax - eps;

// Note that # of "transitions" = # of rows - 1
	double dx = (grid_data.xmax - grid_data.xmin) / (grid_data.grid.size() - 1);
	double dy = (grid_data.ymax - grid_data.ymin) / (grid_data.grid[0].size() - 1);

// indexes of the "square", where (x,y) is located

	int xi = int((x - grid_data.xmin) / dx);
	int yi = int((y - grid_data.ymin) / dy);

// coordinates of (x,y) on the "unit square"
	double xc = (x - grid_data.xmin) / dx - xi;
	double yc = (y - grid_data.ymin) / dy - yi;

//std::printf("%.6f, %.6f, %.6f, %d\n", x, xc, dx, xi);
//std::printf("%.6f, %.6f, %.6f, %d\n", y, yc, dy, yi);

	return grid_data.grid[xi][yi] * (1 - xc) * (1 - yc) + grid_data.grid[xi + 1][yi] * xc * (1 - yc)
			+ grid_data.grid[xi][yi + 1] * (1 - xc) * yc + grid_data.grid[xi + 1][yi + 1] * xc * yc;
}

