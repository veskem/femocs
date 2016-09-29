/*
 * main.cc
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan
 */

#include <cstdlib>
#include <iostream>

#include "laplace.h"
#include "mesh_preparer.h"
#include "physical_quantities.h"
#include "currents_and_heating.h"

#include "field_currents_heating.h"

int main() {


	PhysicalQuantities pq;
	if (!pq.load_emission_data("../res/physical_quantities/gtf_grid_1000x1000.dat"))
			return EXIT_FAILURE;
	if (!pq.load_nottingham_data("../res/physical_quantities/nottingham_grid_300x300.dat"))
			return EXIT_FAILURE;
	if (!pq.load_resistivity_data("../res/physical_quantities/cu_res_mod.dat"))
			return EXIT_FAILURE;



/*
	MeshPreparer<3> mesh_preparer;
	field_currents_heating::FieldCurrentsHeating<3> fch(pq);

	Triangulation<3> *p_mesh = fch.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_mesh, "../res/3d_meshes/mushroom_mod.msh");
	mesh_preparer.mark_top_and_bottom_boundary(p_mesh);
	mesh_preparer.output_mesh(p_mesh, "mesh_3d.vtk");

	fch.run();
*/


	MeshPreparer<2> mesh_preparer;

	laplace::Laplace<2> field;
	Triangulation<2> *p_vmesh = field.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_vmesh, "../res/2d_meshes/vacuum_aligned.msh");
	mesh_preparer.output_mesh(p_vmesh, "vacuum_mesh.vtk");
	field.run();

	currents_heating::CurrentsAndHeating<2> ch(pq, &field);
	Triangulation<2> *p_cmesh = ch.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_cmesh, "../res/2d_meshes/copper_aligned.msh");
	mesh_preparer.output_mesh(p_cmesh, "copper_mesh.vtk");

	ch.run();

	return EXIT_SUCCESS;
}
