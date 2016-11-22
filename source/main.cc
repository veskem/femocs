/*
 * main.cc
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan
 */

#include <cstdlib>
#include <iostream>

#include "laplace.h"
#include "physical_quantities.h"
#include "currents_and_heating.h"


int main() {

	dealii::Timer timer;

	fch::PhysicalQuantities pq;
	if (!pq.load_emission_data("../res/physical_quantities/gtf_grid_1000x1000.dat"))
			return EXIT_FAILURE;
	if (!pq.load_nottingham_data("../res/physical_quantities/nottingham_grid_300x300.dat"))
			return EXIT_FAILURE;
	if (!pq.load_resistivity_data("../res/physical_quantities/cu_res_mod.dat"))
			return EXIT_FAILURE;

	std::cout << "    Load PhysicalQuantities: " << timer.wall_time() << " s" << std::endl; timer.restart();

	/* ------------------------------------------------------------------------------------- */
	/* Laplace solver */
	fch::Laplace<3> laplace_solver;
	laplace_solver.import_mesh_from_file("../res/3d_meshes/nanotip_vacuum.msh");
	laplace_solver.setup_system();
	laplace_solver.assemble_system();
	laplace_solver.solve();
	laplace_solver.output_results("field_sol.vtk");

	std::cout << "    Solved laplace_solver: " << timer.wall_time() << " s" << std::endl; timer.restart();

	// Accessing data in vertexes
	std::vector<int> cell_indexes, vertex_indexes;
	for (int i = 0; i < 10; i++) {
		cell_indexes.push_back(i);
		vertex_indexes.push_back(0);
	}
	std::vector<double> potentials = laplace_solver.get_potential(cell_indexes, vertex_indexes);
	std::vector<dealii::Tensor<1, 3>> fields = laplace_solver.get_field(cell_indexes, vertex_indexes);

	/* ------------------------------------------------------------------------------------- */
	/* Currents and heating */
	fch::CurrentsAndHeating<3> ch_solver(pq, &laplace_solver);
	ch_solver.import_mesh_from_file("../res/3d_meshes/nanotip_copper.msh");

	double final_error = ch_solver.run_specific(10.0, 3, true, "sol", true);

	std::cout << "    Solved currents&heating: " << timer.wall_time() << " s" << std::endl; timer.restart();
	std::cout << "    Final temp. error: " << final_error << std::endl;

	/* ------------------------------------------------------------------------------------- */
	/* Currents and heating resume (or next iteration) */
	fch::CurrentsAndHeating<3> ch_solver2(pq, &laplace_solver, &ch_solver);
	ch_solver2.import_mesh_from_file("../res/3d_meshes/nanotip_copper.msh");

	double final_error2 = ch_solver2.run_specific(10.0, 3, true, "sol_next", true);

	std::cout << "    Solved currents&heating: " << timer.wall_time() << " s" << std::endl; timer.restart();
	std::cout << "    Final temp. error: " << final_error2 << std::endl;

	/* Access data in vertexes */
	std::vector<double> temperatures = ch_solver2.get_temperature(cell_indexes, vertex_indexes);
	std::vector<dealii::Tensor<1, 3>> currents = ch_solver2.get_current(cell_indexes, vertex_indexes);

	for (unsigned int i = 0; i < temperatures.size(); i++) {
		std::cout << i << " T: " << temperatures[i] << "; C: " << currents[i] << std::endl;
	}


/* 2d case usage */
/*
	laplace::Laplace<2> field;
	field.import_mesh_from_file("../res/2d_meshes/vacuum_aligned_dense.msh", "vacuum_mesh_2d.vtk");
	field.run();

	currents_heating::CurrentsAndHeating<2> ch(pq, &field);
	ch.import_mesh_from_file("../res/2d_meshes/copper_aligned_dense.msh", "copper_mesh_2d.vtk");
	ch.run();
*/

/* 2d Mesh splitting */
/*
	MeshPreparer<2> mesh_preparer;
	Triangulation<2> mesh;
	mesh_preparer.import_mesh_from_file(&mesh, "../res/2d_meshes/both_dense.msh");
	Triangulation<2> new_mesh = mesh_preparer.remove_cells_with_id(&mesh, 20);
	mesh_preparer.mark_vacuum_boundary(&new_mesh);
	mesh_preparer.output_mesh(&new_mesh, "vacuum_aligned_dense.msh");
*/

/* FCH usage */
/*
	MeshPreparer<2> mesh_preparer_fch;
	field_currents_heating::FieldCurrentsHeating<2> fch(pq);

	Triangulation<2> *p_mesh = fch.getp_triangulation();
	mesh_preparer_fch.import_mesh_from_file(p_mesh, "../res/2d_meshes/mod.msh");
	mesh_preparer_fch.mark_top_and_bottom_boundary(p_mesh);
	mesh_preparer_fch.output_mesh(p_mesh, "mesh_2d.vtk");

	fch.run();
*/
	return EXIT_SUCCESS;
}



