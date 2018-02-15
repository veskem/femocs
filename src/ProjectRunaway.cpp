/*
 * ProjectRunaway.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include "GeneralProject.h"
#include "ProjectRunaway.h"
#include "Coarseners.h"
#include "Macros.h"
#include "Tethex.h"
#include "VoronoiMesh.h"
#include "Laplace.h"

#include <omp.h>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;
namespace femocs {

ProjectRunaway::ProjectRunaway(AtomReader &a, Config &c) :
        ProjectNanotip(a, c),
        surface_fields(&vacuum_interpolator), surface_temperatures(&bulk_interpolator),
        emission(surface_fields, surface_temperatures, &vacuum_interpolator),
        ch_solver(&ch_solver1), prev_ch_solver(NULL),
        pic_solver(laplace_solver, ch_transient_solver, vacuum_interpolator, emission)
{
    forces.set_interpolator(&vacuum_interpolator);
    temperatures.set_interpolator(&bulk_interpolator);

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    phys_quantities.initialize_with_hc_data();
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_transient_solver.set_physical_quantities(&phys_quantities);
    end_msg(t0);
}

// Write all the available data to file for debugging purposes
int ProjectRunaway::force_output() {
    if (conf.behaviour.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    mesh->hexahedra.write("out/hexmesh_err.vtk");
    mesh->elems.write("out/tetmesh_err.vtk");
    mesh->faces.write("out/trimesh_err.vtk");

    vacuum_interpolator.nodes.write("out/result_E_phi_err.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi_vacuum_err.vtk");

    if (bulk_interpolator.nodes.size() > 0) {
        if (conf.heating.mode == "transient" || conf.heating.mode == "stationary")
            bulk_interpolator.nodes.write("out/result_J_T_err.xyz");
        if (conf.heating.mode == "transient") {
            ch_transient_solver.output_results_current("out/result_J_err.vtk");
            ch_transient_solver.output_results_heating("out/result_T_err.vtk");
        }
        if (conf.heating.mode == "stationary")
            ch_solver->output_results("out/result_J_T_err.vtk");
    }

    return 0;
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int ProjectRunaway::run(const double elfield, const int tstep) {
    stringstream stream;
    stream << fixed << setprecision(3);

    double tstart = omp_get_wtime();

    stream.str("");
    stream << "Atoms haven't moved significantly, " << reader.rms_distance
            << " < " << conf.tolerance.distance << "! Previous mesh will be used!";

    //******************** MESHING *****************************
    if (reinit(tstep)) { // reinit and check skip_meshing
        write_verbose_msg(stream.str());
    } else {
        generate_mesh();
        prepare_fem();
    }

    //****************** RUNNING Field - PIC calculation ********
    skip_meshing = true;

    if (conf.field.solver == "poisson") {
        if (solve_pic(elfield, max(delta_t_MD * 1.e15, conf.pic.total_time))) {
            force_output();
            check_return(true, "Solving PIC failed!");
        }
    } else {
        if (solve_laplace(elfield)) {
            force_output();
            check_return(true, "Solving Laplace equation failed!");
        }
    }

    if (solve_heat(conf.heating.t_ambient)) {
        force_output();
        check_return(true, "Solving heat & continuity equation failed!");
    }

    finalize();

    stream.str("");
    stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());

    return 0;
}

int ProjectRunaway::prepare_fem(){
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(conf.field.E0, conf.tolerance.field_min, conf.tolerance.field_max,
            conf.geometry.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !laplace_solver.import_mesh_directly(mesh->nodes.export_dealii(),
            mesh->hexahedra.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Importing mesh to transient J & T solver...");
    fail = !ch_transient_solver.import_mesh_directly(mesh->nodes.export_dealii(),
            mesh->hexahedra.export_bulk());
    check_return(fail, "Importing bulk mesh to Deal.II failed!");
    end_msg(t0);

    vacuum_interpolator.initialize(mesh);
    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    vacuum_interpolator.lintets.narrow_search_to(TYPES.VACUUM);
    bulk_interpolator.lintets.narrow_search_to(TYPES.BULK);

    return 0;
}

// Run Pic simulation for dt_main time advance
int ProjectRunaway::solve_pic(const double E0, const double dt_main) {

    int time_subcycle = ceil(dt_main / conf.pic.dt_max); // dt_main = delta_t_MD converted to [fs]
    double dt_pic = dt_main/time_subcycle;

    start_msg(t0, "=== Running PIC...\n");
    pic_solver.set_params(conf.field, conf.pic, dt_pic, mesh->nodes.stat);

    // Obtain quadrangle centroids where field will be probed
    surface_fields.set_preferences(false, 2, 3);
    surface_fields.transfer_elfield(ch_transient_solver);

    // Interpolate J & T on quadrangle centroids
    // Temperature do not need to be recalculated, as thermal timestep is >> PIC timestep
    surface_temperatures.set_preferences(false, 2, 3);
    surface_temperatures.interpolate(ch_transient_solver);

    emission.initialize(mesh);

    int n_lost, n_injected, n_cg_steps;
    for (int i = 0; i < time_subcycle; i++) {
        n_lost = pic_solver.update_positions();
        n_cg_steps = pic_solver.run_cycle(i == 0);

        vacuum_interpolator.extract_solution(&laplace_solver);
        surface_fields.calc_interpolation();
        emission.calc_emission(conf.emission, conf.field.V0);

        n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);

        write_verbose_msg("# CG steps|injected|deleted electrons: "
                + to_string(n_cg_steps) + "|" + to_string(n_injected) + "|" + to_string(n_lost));
    }
    end_msg(t0);

    pic_solver.write("out/electrons.movie", 0);

//    const double norm_avg_field = 2.518780;
//    const double norm_max_field = 302.992500;
//    double avg_field = 0;
//    double max_field = 0;
//    for (Solution s : *vacuum_interpolator.nodes.get_solutions()) {
//        avg_field += s.norm;
//        max_field = max(max_field, s.norm);
//    }
//    avg_field /= vacuum_interpolator.nodes.size();
//
//    write_verbose_msg("avg: calc/norm = " + to_string(avg_field) + "/" + to_string(norm_avg_field) + " = " + to_string(avg_field/norm_avg_field));
//    write_verbose_msg("enhancement: calc/norm = " + to_string(max_field/E0) + "/" + to_string(norm_max_field/E0) + " = " + to_string(max_field/norm_max_field));

    return 0;

    //7. Save ions and neutrals that are inbound on the MD domain
    //    somewhere where the MD can find them
    // TODO LATER

    //8. Give the heat- and current fluxes to the temperature solver.
    // TODO LATER
}

// Pick a method to solve heat & continuity equations
int ProjectRunaway::solve_heat(const double T_ambient) {
    if (conf.heating.mode == "stationary") {
        return solve_stationary_heat();
    }
    else if(conf.heating.mode == "transient") {
        double delta_time = delta_t_MD * (timestep - last_full_timestep); //in sec
        int success = solve_transient_heat(delta_time);
        write_verbose_msg("Current and heat advanced for " + to_string(delta_time));
        return success;
    }
    else if (conf.heating.mode == "converge") {
        return solve_converge_heat();
    }

    return 0;
}

// Solve steady-state heat and continuity equations
int ProjectRunaway::solve_stationary_heat() {
    start_msg(t0, "=== Initializing stationary J & T solver...");
    ch_solver->reinitialize(&laplace_solver, prev_ch_solver);
    end_msg(t0);

    stringstream ss; ss << *(ch_solver);
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Importing mesh to J & T solver...");
    fail = !ch_solver->import_mesh_directly(mesh->nodes.export_dealii(),
            mesh->hexahedra.export_bulk());
    check_return(fail, "Importing mesh to Deal.II failed!");
    ch_solver->setup_system();
    end_msg(t0);

    start_msg(t0, "=== Transfering elfield to J & T solver...");
    FieldReader field_reader(&vacuum_interpolator);
    field_reader.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    field_reader.transfer_elfield(ch_solver);
    end_msg(t0);

    field_reader.write("out/surface_field.xyz");

    start_msg(t0, "=== Running J & T solver...\n");
    ch_solver->set_ambient_temperature(conf.heating.t_ambient);
    double t_error = ch_solver->run_specific(conf.heating.t_error, conf.heating.n_newton, false, "", MODES.VERBOSE, 2, 400, true);
    end_msg(t0);

    ch_solver->output_results("out/result_J_T.vtk");

    check_return(t_error > conf.heating.t_error, "Temperature didn't converge, err=" + to_string(t_error));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator.nodes.write("out/result_J_T.xyz");

    // Swap current-and-heat-solvers to use solution from current run as a guess in the next one
    static bool odd_run = true;
    if (odd_run) {
        ch_solver = &ch_solver2; prev_ch_solver = &ch_solver1;
    }
    else {
        ch_solver = &ch_solver1; prev_ch_solver = &ch_solver2;
    }
    odd_run = !odd_run;

    return 0;
}

// Solve transient heat and continuity equations
int ProjectRunaway::solve_transient_heat(const double delta_time) {
    double multiplier = 1.;

    // Transfer elfield to J & T solver
    surface_fields.set_preferences(false, 2, 3);
    surface_fields.transfer_elfield(ch_transient_solver);

    // Interpolate J & T on face centroids
    surface_temperatures.set_preferences(false, 2, 3);
    surface_temperatures.interpolate(ch_transient_solver);

    // Calculate field emission
    emission.initialize(mesh);
    emission.calc_emission(conf.emission, conf.field.V0);
    emission.export_emission(ch_transient_solver);

    start_msg(t0, "=== Setup transient J & T solver...");
    ch_transient_solver.setup_current_system();
    ch_transient_solver.setup_heating_system();
    end_msg(t0);

    start_msg(t0, "=== Calculating current density...");
    ch_transient_solver.assemble_current_system(); // assemble matrix for current density equation; current == electric current
    unsigned int ccg = ch_transient_solver.solve_current();  // ccg == number of current calculation (CG) iterations
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(ccg));

    start_msg(t0, "=== Calculating temperature distribution...");
    ch_transient_solver.set_timestep(delta_time);
    ch_transient_solver.assemble_heating_system_euler_implicit();
    unsigned int hcg = ch_transient_solver.solve_heat(); // hcg == number of temperature calculation (CG) iterations
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(hcg));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    bulk_interpolator.extract_solution(ch_transient_solver);
    end_msg(t0);

    bulk_interpolator.nodes.write("out/result_rho_T.xyz");
    bulk_interpolator.lintets.write("out/result_rho_T.vtk");

    return 0;
}

// Solve transient heat and continuity until convergence is achieved
int ProjectRunaway::solve_converge_heat() {
    static bool first_call = true;

    start_msg(t0, "=== Setup transient J & T solver...");
    emission.initialize(mesh);
    ch_transient_solver.setup_current_system();
    ch_transient_solver.setup_heating_system();
    end_msg(t0);


    double current_time = 0.;
    double delta_time = 1.e-12; //in seconds!!
    double multiplier = 1.;

    for (int i = 0; i < 1000; ++i){

        start_msg(t0, "=== Interpolating J & T on face centroids...");
        surface_temperatures.set_preferences(false, 2, conf.behaviour.interpolation_rank);
        surface_temperatures.interpolate(ch_transient_solver);
        end_msg(t0);
        if (MODES.VERBOSE) surface_temperatures.write("out/surface_temperature.xyz");

        start_msg(t0, "=== Calculating field emission...");
        emission.set_multiplier(multiplier);
        emission.calc_emission(conf.emission, conf.field.V0);
        multiplier = emission.get_multiplier();
        emission.export_emission(ch_transient_solver);
        end_msg(t0);
        if (MODES.VERBOSE) emission.write("out/surface_emission.xyz");

        start_msg(t0, "=== Calculating current density...");
        ch_transient_solver.assemble_current_system(); // assemble matrix for electric current density equation
        unsigned int ccg = ch_transient_solver.solve_current();  // ccg == number of current calculation (CG) iterations
        end_msg(t0);

        start_msg(t0, "=== Calculating temperature distribution...\n");
        ch_transient_solver.set_timestep(delta_time);
        ch_transient_solver.assemble_heating_system_euler_implicit();
        unsigned int hcg = ch_transient_solver.solve_heat(); // hcg == number of temperature calculation (CG) iterations
        end_msg(t0);

        start_msg(t0, "=== Extracting J & T...");
        bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
        bulk_interpolator.extract_solution(ch_transient_solver);
        end_msg(t0);
        bulk_interpolator.nodes.write("out/result_J_T.movie");

        current_time += delta_time;
        if (MODES.VERBOSE) {
            double max_T = ch_transient_solver.get_max_temperature();
            printf("  i=%d ; dt= %5.3e ps ; t= %5.3e ps ; ccg= %2d ; hcg= %2d ; Tmax= %6.2f \n",
                    i, delta_time * 1.e12, current_time * 1.e12, ccg, hcg, max_T);
        }

        if (max(hcg, ccg) < 120 || hcg < 30)
            delta_time *= 1.25;
        else if (max(hcg, ccg) > 150)
            delta_time /= 1.25;

        if (max(hcg, ccg) < 10) return 0;
    }
    write_silent_msg("WARNING: Heat equation did not converged after 1000 steps.");

    return 0;
}

} /* namespace femocs */
