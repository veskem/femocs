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

	MeshPreparer<2> mesh_preparer;
	PhysicalQuantities pq;
	if (!pq.load_emission_data("../res/physical_quantities/gtf_grid_1000x1000.dat"))
			return EXIT_FAILURE;
	if (!pq.load_nottingham_data("../res/physical_quantities/nottingham_grid_300x300.dat"))
			return EXIT_FAILURE;
	if (!pq.load_resistivity_data("../res/physical_quantities/cu_res_mod.dat"))
			return EXIT_FAILURE;
	field_currents_heating::FieldCurrentsHeating<2> fch(pq);
	Triangulation<2> *p_mesh = fch.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_mesh, "../res/2d_meshes/both.msh");
	mesh_preparer.output_mesh(p_mesh, "both_mesh_2d.vtk");

	fch.run();

/*
	MeshPreparer<2> mesh_preparer;
	Laplace<2> laplace_problem;
	Triangulation<2> *p_vacuum_mesh = laplace_problem.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_vacuum_mesh, "../res/2d_meshes/both.msh");
	mesh_preparer.output_mesh(p_vacuum_mesh, "both_mesh_2d.vtk");
*/

/*
	Timer timer;
	timer.start();
	PhysicalQuantities pq;
	if (!pq.load_emission_data("../res/physical_quantities/gtf_grid_1000x1000.dat"))
		return EXIT_FAILURE;
	if (!pq.load_resistivity_data("../res/physical_quantities/cu_res_mod.dat"))
		return EXIT_FAILURE;
	timer.stop();
	std::cout << "-------- Loaded physical quantities: " << timer() << "s" << std::endl;

	std::cout << pq.evaluate_current(1.0, 1000.0) << std::endl;
*/
/*
	MeshPreparer<3> mesh_preparer;
	Laplace<3> laplace_problem;
	Timer timer;

	timer.start ();
	Triangulation<3> *p_vacuum_mesh = laplace_problem.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_vacuum_mesh, "../res/tethex_vacuum.msh");
	mesh_preparer.output_mesh(p_vacuum_mesh, "vacuum_mesh.vtk");
	mesh_preparer.mark_vacuum_boundary(p_vacuum_mesh);
	timer.stop();
	std::cout << "-------- Imported, outputted and marked vacuum mesh: " << timer() << "s" << std::endl;


	timer.restart();
	laplace_problem.run();
	timer.stop();
	std::cout << "-------- Solved Laplace problem: " << timer() << "s" << std::endl;
*/
/*
	timer.restart();
	for (int i = 0; i < 100; i++) {
		laplace_problem.probe_field(Point<3>(20.0, 20.0, 40.0));
	}
	timer.stop();
	std::cout << "Probing field: " << timer() << std::endl;
*/

/*
	timer.restart();
	PhysicalQuantities physical_quantities;
	if (!physical_quantities.load_emission_data("../res/physical_quantities/gtf_grid_1000x1000.dat"))
		return EXIT_FAILURE;
	if (!physical_quantities.load_resistivity_data("../res/physical_quantities/cu_res_mod.dat"))
		return EXIT_FAILURE;
	timer.stop();
	std::cout << "-------- Loaded physical quantities: " << timer() << "s" << std::endl;

	timer.restart();
	CurrentsAndHeating<3> currents_and_heating(physical_quantities, &laplace_problem);
	Triangulation<3> *p_copper_mesh = currents_and_heating.getp_triangulation();
	mesh_preparer.import_mesh_from_file(p_copper_mesh, "../res/tethex_copper.msh");
	mesh_preparer.output_mesh(p_copper_mesh, "copper_mesh.vtk");
	mesh_preparer.mark_copper_boundary(p_copper_mesh);
	timer.stop();
	std::cout << "-------- Imported, outputted and marked copper mesh: " << timer() << "s" << std::endl;

	timer.restart();
	currents_and_heating.run();
	timer.stop();
	std::cout << "-------- Solved currents and heating: " << timer() << "s" << std::endl;
*/

	return EXIT_SUCCESS;
}
