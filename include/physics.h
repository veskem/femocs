/*
 * physics.h
 *
 *  Created on: Jan 15, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_PHYSICS_H_
#define INCLUDE_PHYSICS_H_

#include <cmath>

namespace Emitter {

double emission_current(double field, double temp) {
	double j_e = field*field*std::exp(-1.0/field)/(80.0*3.0*3.0); // field emission
	double j_t = temp*temp/(1000.0*1500.0*1500.0); //thermionic emission
	return j_e+j_t;
}

/**
 * Electrical conductivity
 */
double el_conductivity(double temperature) {
	return 5/std::max(temperature,300.0);
}

/**
 * Thermal conductivity
 */
double th_conductivity(double temperature) {
	return 0.5/std::max(temperature,300.0);
}

}

#endif /* INCLUDE_PHYSICS_H_ */
