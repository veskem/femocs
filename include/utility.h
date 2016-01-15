/*
 * utility.h
 *
 *  Created on: Jan 15, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_UTILITY_H_
#define INCLUDE_UTILITY_H_


#include <deal.II/grid/tria.h>

#include <string>

namespace Emitter {

using namespace dealii;

// Method to output the mesh to a specified .eps file
void output_mesh(Triangulation<2>& mesh, std::string name);

}

#endif /* INCLUDE_UTILITY_H_ */
