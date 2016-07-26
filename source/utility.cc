/*
 * utility.cc
 *
 *  Created on: Jan 15, 2016
 *      Author: kristjan
 */


#include <deal.II/grid/grid_out.h>

#include <fstream>

#include "utility.h"


void output_mesh(Triangulation<2>& mesh, std::string name) {
	std::ofstream out(name);
	GridOut grid_out;
	grid_out.write_eps(mesh, out);
	std::cout << "Grid written to " << name << std::endl;
}



