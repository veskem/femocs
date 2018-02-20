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
#include "DealSolver.h"

#include <omp.h>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;
namespace femocs {

ProjectRunaway::ProjectRunaway(AtomReader &a, Config &conf_) :
        ProjectNanotip(a, conf_),
        bulk_interpolator("rho", "temperature"),
        surface_fields(&vacuum_interpolator),
        surface_temperatures(&bulk_interpolator),
        emission(&surface_fields, &surface_temperatures, &vacuum_interpolator),
        phys_quantities(conf.heating),
        pic_solver(&poisson_solver, &ch_solver, &vacuum_interpolator, &emission)
{
    forces.set_interpolator(&vacuum_interpolator);
    temperatures.set_interpolator(&bulk_interpolator);
    poisson_solver.set_particles(pic_solver.get_particles());

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    phys_quantities.initialize_with_hc_data();
    ch_solver.set_dependencies(&phys_quantities, &conf_.heating);
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
        if (conf.heating.mode == "transient") {
            bulk_interpolator.nodes.write("out/result_J_T_err.xyz");
            ch_solver.current.write("out/result_J_err.vtk");
            ch_solver.heat.write("out/result_T_err.vtk");
        }
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
    }

    //****************** RUNNING Field - PIC calculation ********
    skip_meshing = true;

    if (conf.field.solver == "poisson") {
        double delta_t = max(delta_t_MD * 1.e15, conf.pic.total_time);
        if (solve_pic(elfield, delta_t)) {
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

// Run Pic simulation for dt_main time advance
int ProjectRunaway::solve_pic(const double E0, const double dt_main) {

    int time_subcycle = ceil(dt_main / conf.pic.dt_max); // dt_main = delta_t_MD converted to [fs]
    double dt_pic = dt_main/time_subcycle;
    
    start_msg(t0, "=== Importing mesh to Poisson solver...");
    fail = !poisson_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexahedra.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);
    
    start_msg(t0, "=== Initializing Poisson solver...");
    poisson_solver.setup(-E0, conf.field.V0);
    end_msg(t0);

    start_msg(t0, "=== Importing mesh to J & T solver...");
    fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexahedra.export_bulk());
    check_return(fail, "Importing bulk mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Initializing PIC dependencies...");
    vacuum_interpolator.initialize(mesh);
    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    vacuum_interpolator.lintets.narrow_search_to(TYPES.VACUUM);
    bulk_interpolator.lintets.narrow_search_to(TYPES.BULK);
    pic_solver.set_params(conf.field, conf.pic, dt_pic, mesh->nodes.stat);
    end_msg(t0);

    start_msg(t0, "=== Interpolating elfield on faces...");
    surface_fields.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    surface_fields.interpolate(ch_solver);
    end_msg(t0);
    
    surface_fields.write("out/surface_field.xyz");
    
    // Temperature do not need to be recalculated, as thermal timestep is >> PIC timestep
    start_msg(t0, "=== Interpolating J & T on faces...");
    surface_temperatures.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    surface_temperatures.interpolate(ch_solver);
    end_msg(t0);
    
    surface_temperatures.write("out/surface_temperature.xyz");

    emission.initialize(mesh);
    
    start_msg(t0, "=== Running PIC...\n");
    int n_lost, n_injected, n_cg_steps;
    
    for (int i = 0; i < time_subcycle; i++) {
        n_lost = pic_solver.update_positions();
        n_cg_steps = pic_solver.run_cycle(i == 0);

        vacuum_interpolator.extract_solution(poisson_solver);
        surface_fields.calc_interpolation();
        emission.calc_emission(conf.emission, conf.field.V0);

        n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);
        
        if (MODES.VERBOSE)
            printf("  #CG steps=%d, max field=%.3f, #injected|deleted electrons=%d|%d\n",
                    n_cg_steps, max_field(), n_injected, n_lost);
    }
    
    end_msg(t0);

    pic_solver.write("out/electrons.movie", 0);

    return 0;

    //7. Save ions and neutrals that are inbound on the MD domain somewhere where the MD can find them
    // TODO LATER
    //8. Give the heat- and current fluxes to the temperature solver.
    // TODO LATER
}

double ProjectRunaway::max_field() {
    double max_field = 0;
    for (Solution s : *vacuum_interpolator.nodes.get_solutions())
        max_field = max(max_field, s.norm);
    return max_field;
}

// Pick a method to solve heat & continuity equations
int ProjectRunaway::solve_heat(const double T_ambient) {
    if(conf.heating.mode == "transient") {
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

// Solve transient heat and continuity equations
int ProjectRunaway::solve_transient_heat(const double delta_time) {
    double multiplier = 1.;

    // Interpolate elfield on face centroids
    surface_fields.set_preferences(false, 2, 3);
    surface_fields.interpolate(ch_solver);

    // Interpolate J & T on face centroids
    surface_temperatures.set_preferences(false, 2, 3);
    surface_temperatures.interpolate(ch_solver);

    // Calculate field emission
    emission.initialize(mesh);
    emission.calc_emission(conf.emission, conf.field.V0);
    emission.export_emission(ch_solver);

    start_msg(t0, "=== Setup transient J & T solver...");
    ch_solver.setup(conf.heating.t_ambient);
    end_msg(t0);

    start_msg(t0, "=== Calculating current density...");
    ch_solver.current.assemble();
    unsigned int ccg = ch_solver.current.solve();
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(ccg));

    start_msg(t0, "=== Calculating temperature distribution...");
    ch_solver.heat.assemble(delta_time);
    unsigned int hcg = ch_solver.heat.solve();
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(hcg));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator.nodes.write("out/result_J_T.xyz");
    bulk_interpolator.lintets.write("out/result_J_T.vtk");

    return 0;
}

// Solve transient heat and continuity until convergence is achieved
int ProjectRunaway::solve_converge_heat() {

    start_msg(t0, "=== Importing mesh to J & T solver...");
    fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexahedra.export_bulk());
    check_return(fail, "Importing bulk mesh to Deal.II failed!");
    end_msg(t0);

    bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
    bulk_interpolator.lintets.narrow_search_to(TYPES.BULK);

    start_msg(t0, "=== Interpolating elfield on faces...");
    surface_fields.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    surface_fields.interpolate(ch_solver);
    end_msg(t0);
    surface_fields.write("out/surface_field.xyz");

    start_msg(t0, "=== Setup transient J & T solver...");
    ch_solver.setup(conf.heating.t_ambient);
    end_msg(t0);

    surface_temperatures.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    emission.initialize(mesh);

    double current_time = 0.;
    double delta_time = 1.e-12; //in seconds!!
    double multiplier = 1.;

    start_msg(t0, "=== Running converge J & T solver...\n");
    int converge_steps = 0;
    for (; converge_steps < 1000; ++converge_steps) {
        double tstart = omp_get_wtime();

        // Interpolate J & T on face centroids
        surface_temperatures.interpolate(ch_solver);

        // Calculating field emission
        multiplier = emission.calc_emission(multiplier, conf.emission, conf.field.V0);
        emission.export_emission(ch_solver);

        // Calculate current density
        ch_solver.current.assemble();
        unsigned int ccg = ch_solver.current.solve(); // ccg == nr of CG iterations

        // Calculate temperature distribution
        ch_solver.heat.assemble(delta_time);
        unsigned int hcg = ch_solver.heat.solve(); // hcg == nr of CG iterations

        // Extract J & T
        bulk_interpolator.extract_solution(ch_solver);

        if (MODES.VERBOSE) {
            double max_T = ch_solver.heat.max_solution();
            double time = omp_get_wtime() - tstart;
            printf("  i=%d, t=%.2fps, dt=%.2fps, ccg=%d, hcg=%d, time=%.3fs, Tmax=%.2fK\n",
                    converge_steps, current_time * 1.e12, delta_time * 1.e12, ccg, hcg, time, max_T);
        }
        current_time += delta_time;

        if (max(hcg, ccg) < 120 || hcg < 30)
            delta_time *= 1.25;
        else if (max(hcg, ccg) > 150)
            delta_time /= 1.25;

        if (max(hcg, ccg) < 10) break;
    }
    end_msg(t0);
    if (converge_steps >= 1000)
        write_silent_msg("WARNING: Heat equation did not converge after 1000 steps!");

    surface_temperatures.write("out/surface_temperature.xyz");
    emission.write("out/surface_emission.xyz");
    bulk_interpolator.nodes.write("out/result_J_T.xyz");

    return 0;
}

} /* namespace femocs */
