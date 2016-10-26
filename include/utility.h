/*
 * utility.h
 *
 *  Created on: Oct 26, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_UTILITY_H_
#define INCLUDE_UTILITY_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

using namespace dealii;

template <int dim>
unsigned nearest_point_index(const Point<dim> p, const std::vector<Point<dim>> points, bool print = false) {
	double best_distance = 1e16;
	unsigned best_i = 0;
	for (unsigned i = 0; i<points.size(); i++) {
		double dist = p.distance(points[i]);
		if (dist < best_distance) {
			best_distance = dist;
			best_i = i;
		}
	}
	if (print) std::cout << best_distance << std::endl;
	return best_i;
}


#endif /* INCLUDE_UTILITY_H_ */
