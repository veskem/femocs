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

    fch::Laplace<2> laplace;
    laplace.import_mesh_from_file("../res/2d_meshes/vacuum_aligned.msh");
    laplace.set_applied_efield(8.0);
    laplace.run();

    double time_step = 0.01e-12; // seconds
    fch::CurrentsAndHeating<2> ch(time_step, &pq);
    ch.import_mesh_from_file("../res/2d_meshes/copper_aligned.msh");

    ch.setup_current_system();
    ch.setup_heating_system();
    ch.set_electric_field_bc(laplace);

    for (int i = 0; i < 1; i++) {
        ch.assemble_current_system();
        unsigned int ccg = ch.solve_current();

        ch.assemble_heating_system();
        unsigned int hcg = ch.solve_heat();
        double max_T = ch.get_max_temperature();
        std::printf("    t=%4.2fps; ccg=%2d; hcg=%2d; max_T=%6.2f\n", i*time_step*1e12, ccg, hcg, max_T);

        if (i%1 == 0) {
            ch.output_results_current("./output/current_solution-"+std::to_string(i)+".vtk");
            ch.output_results_heating("./output/heat_solution-"+std::to_string(i)+".vtk");
        }
    }

    return EXIT_SUCCESS;
}

