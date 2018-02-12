/*
 * main.cc -> Test_solver.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan, Veske
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/stat.h>  // for checking if directory exists

#include <deal.II/base/timer.h>

#include "CurrentsAndHeating.h"
#include "CurrentsAndHeatingStationary.h"
#include "Laplace.h"
#include "PhysicalQuantities.h"

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
        std::cout
                << "res/ folder not found in " + res_path1 + " or " + res_path2
                        + ". Exiting..." << std::endl;
        return EXIT_FAILURE;
    }

    // If physical quantities are not read from a file, hardcoded values will be used.
    if (!(pq.load_emission_data(
            res_path + "/physical_quantities/gtf_200x200.dat")
            && pq.load_nottingham_data(
                    res_path + "/physical_quantities/nottingham_200x200.dat")
            && pq.load_resistivity_data(
                    res_path + "/physical_quantities/cu_res.dat"))) {
        std::cout << "Couldn't load pq data, using default values..."
                << std::endl;
    } else {
        std::cout << "    Loaded PhysicalQuantities: " << timer.wall_time()
                << " s" << std::endl;
        timer.restart();
    }


// Transient example

    fch::Laplace<2> laplace;
    laplace.import_mesh_from_file("../res/2d_meshes/vacuum_aligned.msh");
    laplace.set_applied_efield(10.0);
    laplace.run();

    double time_step = 0.1e-15; // seconds
    fch::CurrentsAndHeating<2> ch(time_step, &pq);
    ch.import_mesh_from_file("../res/2d_meshes/copper_aligned.msh");

    ch.setup_current_system();
    ch.setup_heating_system();
    //ch.set_electric_field_bc(13.0);
    ch.set_electric_field_bc(laplace);

    int i = 0;
    for (double time = 0.0; time <= 3.0e-15; ) {
        time+=time_step;

        ch.assemble_current_system();
        unsigned int ccg = ch.solve_current();

        if (i == 0) {
            ch.assemble_heating_system_euler_implicit();
        } else {
            //ch.assemble_heating_system_euler_implicit();
            ch.assemble_heating_system_crank_nicolson();
        }

        unsigned int hcg = ch.solve_heat();
        double max_T = ch.get_max_temperature();
        std::printf("    t=%5.3ffs; ccg=%2d; hcg=%2d; max_T=%6.2f\n", time*1e15, ccg, hcg, max_T);

        if (i%10 == 0) {
            ch.output_results_current("./output/current_solution-"+std::to_string(i)+".vtk");
            ch.output_results_heating("./output/heat_solution-"+std::to_string(i)+".vtk");
        }
        i++;
    }


// Simple Stationary 3d usage //
/*
    fch::Laplace<3> laplace_solver;
    laplace_solver.import_mesh_from_file(res_path + "/3d_meshes/vacuum_0.msh");
    laplace_solver.set_applied_efield(1.5);
    laplace_solver.run();

    fch::CurrentsAndHeatingStationary<3> ch_solver;
    ch_solver.reinitialize(&laplace_solver);
    ch_solver.set_physical_quantities(&pq);
    ch_solver.import_mesh_from_file(res_path + "/3d_meshes/copper_0.msh");

    ch_solver.setup_system();
    ch_solver.run_specific(1.0, 100, true, "output/sol", true, 2.0);
*/

// Stationary 3d Test usage with interpolation //
/*
    fch::CurrentsAndHeatingStationary<3> ch_solver1;
    fch::CurrentsAndHeatingStationary<3> ch_solver2;
    fch::CurrentsAndHeatingStationary<3>* ch_solver = &ch_solver1;
    fch::CurrentsAndHeatingStationary<3>* prev_ch_solver = NULL;

    bool even = true;

    ch_solver1.set_physical_quantities(&pq);
    ch_solver2.set_physical_quantities(&pq);

    for (int n = 0; n <= 3; n++) {

        std::cout << "Iteration " << n << std::endl;

        fch::Laplace<3> laplace_solver;
        laplace_solver.set_applied_efield(1.5);
        laplace_solver.import_mesh_from_file(
                res_path + "/3d_meshes/vacuum_" + std::to_string(n) + ".msh");
        double mesh_imp_time = timer.wall_time();
        timer.restart();
        laplace_solver.output_mesh(
                "output/vacuum_" + std::to_string(n) + ".vtk");
        double mesh_exp_time = timer.wall_time();
        timer.restart();
        laplace_solver.setup_system();
        double setup_time = timer.wall_time();
        timer.restart();
        laplace_solver.assemble_system();
        double assemble_time = timer.wall_time();
        timer.restart();
        laplace_solver.solve(2000, 1e-9, true, 1.2);
        double solution_time = timer.wall_time();
        timer.restart();
        laplace_solver.output_results(
                "output/field_sol_" + std::to_string(n) + ".vtk");
        double outp_time = timer.wall_time();
        timer.restart();

        std::cout << "    laplace_solver info:\n        " << laplace_solver
                << std::endl;
        printf(
                "    Timings: mesh_imp: %5.2f; mesh_exp: %5.2f; setup: %5.2f; assemble: %5.2f; sol: %5.2f; outp: %5.2f\n",
                mesh_imp_time, mesh_exp_time, setup_time, assemble_time,
                solution_time, outp_time);

        ch_solver->reinitialize(&laplace_solver, prev_ch_solver); // with IC interpolation
        //ch_solver->reinitialize(&laplace_solver); // without IC interpolation
        ch_solver->import_mesh_from_file(
                res_path + "/3d_meshes/copper_" + std::to_string(n) + ".msh");
        double c_mesh_imp_time = timer.wall_time();
        timer.restart();
        ch_solver->output_mesh("output/copper_" + std::to_string(n) + ".vtk");
        double c_mesh_exp_time = timer.wall_time();
        timer.restart();

        ch_solver->setup_system();

        std::cout << "    ch_solver info:\n        " << *ch_solver << std::endl;
        printf("    Timings: mesh_imp: %5.2f; mesh_exp: %5.2f\n",
                c_mesh_imp_time, c_mesh_exp_time);

        double final_error = ch_solver->run_specific(1.0, 100, true,
                "output/sol_" + std::to_string(n), true, 2.0);

        std::cout << "    Solved currents&heating: " << timer.wall_time()
                << " s" << std::endl;
        timer.restart();
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
*/

// 3d mushroom //
    /*
     fch::Laplace<3> laplace_solver;
     laplace_solver.set_applied_efield(1.0);
     laplace_solver.import_mesh_from_file(res_path+"/3d_meshes/mushroom_vacuum.msh");
     laplace_solver.output_mesh("output/mushroom_vacuum.msh");
     laplace_solver.setup_system();
     laplace_solver.assemble_system();
     laplace_solver.solve();
     laplace_solver.output_results("output/field_sol.vtk");

     fch::CurrentsAndHeatingStationary<3> ch_solver(&pq, &laplace_solver);
     ch_solver.import_mesh_from_file(res_path+"/3d_meshes/mushroom_copper.msh");
     ch_solver.setup_system();

     ch_solver.run_specific(1.0, 100, true, "output/sol", true, 2.0);
     */

// 2d case usage //
    /*
     fch::Laplace<2> laplace;
     laplace.import_mesh_from_file("../res/2d_meshes/vacuum_aligned.msh");
     laplace.output_mesh("output/vacuum_mesh_2d.vtk");
     laplace.set_applied_efield(12.0);
     laplace.run();

     fch::CurrentsAndHeatingStationary<2> ch(&pq, &laplace);
     ch.import_mesh_from_file("../res/2d_meshes/copper_aligned.msh");
     ch.output_mesh("output/copper_mesh_2d.vtk");
     ch.run();
     */

// 2d Mesh splitting //
    /*
     fch::MeshPreparer<2> mesh_preparer;
     dealii::Triangulation<2> mesh;
     mesh_preparer.import_mesh_from_file(&mesh, "../res/2d_meshes/simple.msh");
     dealii::Triangulation<2> new_mesh = mesh_preparer.remove_cells_with_id(&mesh, 1);
     mesh_preparer.mark_copper_boundary(&new_mesh);
     mesh_preparer.output_mesh(&new_mesh, "simple_copper.msh");
     */

// Merged mesh fch usage //
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

