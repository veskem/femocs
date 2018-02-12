/*
 * utility.h
 *
 *  Created on: Oct 26, 2016
 *      Author: kristjan
 *
 *  Some useful snippet tools
 */

#ifndef INCLUDE_UTILITY_H_
#define INCLUDE_UTILITY_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>

namespace fch {

using namespace dealii;

template<int dim>
unsigned nearest_point_index(const Point<dim> p,
        const std::vector<Point<dim>> points, bool print = false) {
    double best_distance = 1e16;
    unsigned best_i = 0;
    for (unsigned i = 0; i < points.size(); i++) {
        double dist = p.distance(points[i]);
        if (dist < best_distance) {
            best_distance = dist;
            best_i = i;
        }
    }
    if (print)
        std::cout << best_distance << std::endl;
    return best_i;
}

template<typename T>
T vector_median(std::vector<T> vector) {
    std::nth_element(vector.begin(), vector.begin() + vector.size() / 2,
            vector.end());
    T upper_median = vector[vector.size() / 2];

    if (vector.size() % 2 == 0) {
        std::nth_element(vector.begin(), vector.begin() + vector.size() / 2 - 1,
                vector.end());
        T lower_median = vector[vector.size() / 2 - 1];
        return (upper_median + lower_median) / 2.0;
    }

    return upper_median;
}

template<typename T>
T vector_mean(const std::vector<T> vector) {
    T sum = 0.0;
    for (T elem : vector) {
        sum += elem;
    }
    return sum / vector.size();
}

template<typename T>
T vector_stdev(const std::vector<T> vector) {
    T mean = vector_mean(vector);
    T sq_sum = 0.0;
    for (T elem : vector) {
        sq_sum += (elem - mean) * (elem - mean);
    }
    return std::sqrt(sq_sum / vector.size());
}

inline bool contains_digit(const std::string& s) {
    std::string::const_iterator it = s.begin();
    for (; it < s.end(); it++) {
        if (std::isdigit(*it))
            return true;
    }
    return false;
}

} // namespace fch

#endif /* INCLUDE_UTILITY_H_ */
