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
	MeshPreparer<2> mesh_preparer_fch;
	field_currents_heating::FieldCurrentsHeating<2> fch(pq);

	Triangulation<2> *p_mesh = fch.getp_triangulation();
	mesh_preparer_fch.import_mesh_from_file(p_mesh, "../res/2d_meshes/mod.msh");
	mesh_preparer_fch.mark_top_and_bottom_boundary(p_mesh);
	mesh_preparer_fch.output_mesh(p_mesh, "mesh_2d.vtk");

	fch.run();
*/

/*
	MeshPreparer<2> mesh_preparer;

	laplace::Laplace<2> field;
	Triangulation<2> *p_vmesh = field.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_vmesh, "../res/2d_meshes/vacuum_aligned_dense.msh");
	mesh_preparer.output_mesh(p_vmesh, "vacuum_mesh.vtk");
	field.run();

	currents_heating::CurrentsAndHeating<2> ch(pq, &field);
	Triangulation<2> *p_cmesh = ch.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_cmesh, "../res/2d_meshes/copper_aligned_dense.msh");
	mesh_preparer.output_mesh(p_cmesh, "copper_mesh.vtk");
	ch.run();
*/

	MeshPreparer<3> mesh_preparer;

	laplace::Laplace<3> field;
	Triangulation<3> *p_vmesh = field.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_vmesh, "../res/3d_meshes/mushroom_vacuum.msh");
	mesh_preparer.mark_vacuum_boundary(p_vmesh);
	mesh_preparer.output_mesh(p_vmesh, "vacuum_mesh.vtk");
	field.run();

	currents_heating::CurrentsAndHeating<3> ch(pq, &field, 898, 0.0137);
	Triangulation<3> *p_cmesh = ch.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_cmesh, "../res/3d_meshes/mushroom_copper.msh");
	mesh_preparer.mark_copper_boundary(p_cmesh);
	mesh_preparer.output_mesh(p_cmesh, "copper_mesh.vtk");
	ch.run();


/* Mesh creation */
/*
	MeshPreparer<2> mesh_preparer;
	Triangulation<2> mesh;
	mesh_preparer.import_mesh_from_file(&mesh, "../res/2d_meshes/both_dense.msh");
	Triangulation<2> new_mesh = mesh_preparer.remove_cells_with_id(&mesh, 20);
	mesh_preparer.mark_vacuum_boundary(&new_mesh);
	mesh_preparer.output_mesh(&new_mesh, "vacuum_aligned_dense.msh");
*/
	return EXIT_SUCCESS;
}



