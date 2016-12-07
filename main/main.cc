/*
 * main.cc
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/stat.h>  // for checking if directory exists

#include <deal.II/base/timer.h>

#include "laplace.h"
#include "physical_quantities.h"
#include "currents_and_heating.h"


int main() {

	dealii::Timer timer;

	fch::PhysicalQuantities pq;

	std::string res_path1 = "../res";
	std::string res_path2 = "heating/res";

	std::string res_path = "";

	struct stat info;
	if (stat(res_path1.c_str(), &info) == 0) {
		res_path = res_path1;
	} else if (stat(res_path2.c_str(), &info) == 0) {
		res_path = res_path2;
	} else {
		std::cout << "res/ folder not found in " + res_path1 + " or " + res_path2 + ". Exiting..." << std::endl;
		return EXIT_FAILURE;
	}

	if (!(pq.load_emission_data(res_path+"/physical_quantities/gtf_grid_1000x1000.dat") &&
		  pq.load_nottingham_data(res_path+"/physical_quantities/nottingham_grid_300x300.dat") &&
		  pq.load_resistivity_data(res_path+"/physical_quantities/cu_res.dat"))) {
		std::cout << "Couldn't load pq data, exiting..." << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "    Loaded PhysicalQuantities: " << timer.wall_time() << " s" << std::endl; timer.restart();

    fch::CurrentsAndHeating<3> ch_solver1;
    fch::CurrentsAndHeating<3> ch_solver2;
    fch::CurrentsAndHeating<3>* ch_solver = &ch_solver1;
    fch::CurrentsAndHeating<3>* prev_ch_solver = NULL;

    bool even = true;

    ch_solver1.set_physical_quantities(&pq);
    ch_solver2.set_physical_quantities(&pq);

	for (int n = 0; n<=3; n++) {

		std::cout << "Iteration " << n << std::endl;

		fch::Laplace<3> laplace_solver;
		laplace_solver.set_applied_efield(1.5);
		laplace_solver.import_mesh_from_file(res_path+"/3d_meshes/vacuum_"+std::to_string(n)+".msh");
		double mesh_imp_time = timer.wall_time(); timer.restart();
		laplace_solver.output_mesh("output/vacuum_"+std::to_string(n)+".vtk");
		double mesh_exp_time = timer.wall_time(); timer.restart();
		laplace_solver.setup_system();
		double setup_time = timer.wall_time(); timer.restart();
		laplace_solver.assemble_system();
		double assemble_time = timer.wall_time(); timer.restart();
		laplace_solver.solve();
		double solution_time = timer.wall_time(); timer.restart();
		laplace_solver.output_results("output/field_sol_"+std::to_string(n)+".vtk");
		double outp_time = timer.wall_time(); timer.restart();

		std::cout << "    laplace_solver info:\n        " << laplace_solver << std::endl;
		printf("    Timings: mesh_imp: %5.2f; mesh_exp: %5.2f; setup: %5.2f; assemble: %5.2f; sol: %5.2f; outp: %5.2f\n",
				mesh_imp_time, mesh_exp_time, setup_time, assemble_time, solution_time, outp_time);
		//std::cout << "    Solved laplace_solver: " << timer.wall_time() << " s" << std::endl; timer.restart();


		ch_solver->reinitialize(&laplace_solver, prev_ch_solver); // with IC interpolation
		//ch_solver->reinitialize(&laplace_solver); // without IC interpolation
		ch_solver->import_mesh_from_file(res_path+"/3d_meshes/copper_"+std::to_string(n)+".msh");
		double c_mesh_imp_time = timer.wall_time(); timer.restart();
		ch_solver->output_mesh("output/copper_"+std::to_string(n)+".vtk");
		double c_mesh_exp_time = timer.wall_time(); timer.restart();

		std::cout << "    Currents and heating:" << std::endl;
		printf("    Timings: mesh_imp: %5.2f; mesh_exp: %5.2f\n", c_mesh_imp_time, c_mesh_exp_time);

		double final_error = ch_solver->run_specific(1.0, 100, true, "output/sol_"+std::to_string(n), true, 1.0);

		std::cout << "    ch_solver info:\n        " << *ch_solver << std::endl;
		std::cout << "    Solved currents&heating: " << timer.wall_time() << " s" << std::endl; timer.restart();
		std::cout << "    Final temp. error: " << final_error << std::endl;

		if (even) {
			ch_solver = &ch_solver2;
			prev_ch_solver = &ch_solver1;
			even = false;
		} else {
			ch_solver = &ch_solver1;
			prev_ch_solver = &ch_solver2;
			even = true;
		}
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



