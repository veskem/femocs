/*
 * mesh_input.h
 *
 *  Created on: Jul 27, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_MESH_PREPARER_H_
#define INCLUDE_MESH_PREPARER_H_


#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>



using namespace std;
using namespace dealii;

enum BoundaryId {
	copper_surface = 0,
	vacuum_top = 1,
	copper_bottom = 2,
	vacuum_sides = 3,
	copper_sides = 4
};

template <int dim>
class MeshPreparer {


	const std::string get_file_ext(const std::string file_name);

public:

	void import_mesh_from_file(Triangulation<dim> *triangulation, const std::string file_name);

	void output_mesh(Triangulation<dim> *triangulation, std::string name);

	void mark_vacuum_boundary(Triangulation<dim> *triangulation);

	void mark_copper_boundary(Triangulation<dim> *triangulation);

	void mark_top_and_bottom_boundary(Triangulation<dim> *triangulation);


};


#endif /* INCLUDE_MESH_PREPARER_H_ */
