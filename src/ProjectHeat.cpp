/*
 * ProjectHeat.cpp
 *
 *  Created on: 4.5.2018
 *      Author: kyritsak
 */

#include "ProjectHeat.h"

namespace femocs {

ProjectHeat::ProjectHeat(AtomReader &reader, Config &conf) : ProjectRunaway(reader, conf)
{}

int ProjectHeat::run(int timestep, double time) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (generate_mesh())
        return process_failed("Mesh generation failed!");

    check_return(!mesh_changed, "First meshing failed! Terminating...");

    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    //***** Run FEM solvers *****

    for(auto factor : conf.SC.apply_factors){
        conf.field.E0 *= factor;
        conf.field.V0 *= factor;

        if (run_field_solver())
            return process_failed("Running field solver in a " + conf.field.mode + " mode failed!");

        if (run_heat_solver())
            return process_failed("Running heat solver in a " + conf.heating.mode + " mode failed!");

        //***** Prepare for data export and next run *****
        if (prepare_export())
            return process_failed("Interpolating solution on atoms failed!");

        finalize(tstart, time);
        conf.field.E0 /= factor;
        conf.field.V0 /= factor;
    }

    return 0;
}

int ProjectHeat::run_field_solver() {
    if (conf.field.mode == "transient")
        return solve_pic(conf.behaviour.timestep_fs, mesh_changed);
    else if (conf.field.mode == "converge")
        return converge_pic(1.e4);
    else if (mesh_changed)
        return solve_laplace(conf.field.E0, conf.field.V0);

    return 0;
}

int ProjectHeat::run_heat_solver() {
    int ccg, hcg;

    if (conf.heating.mode == "converge")
        return converge_heat(conf.heating.t_ambient);

    if (mesh_changed && conf.heating.mode == "transient")
        return solve_heat(conf.heating.t_ambient, GLOBALS.TIME - last_heat_time, true, ccg, hcg);

    return 0;
}

int ProjectHeat::converge_pic(double max_time) {
    double time_window; //time window to check convergence
    int i_max; //window iterations
    if (max_time < conf.pic.dt_max * 16) {
        time_window = max_time;
        i_max = 1;
    } else {
        i_max =  ceil(max_time / (16 * conf.pic.dt_max));
        time_window = max_time / i_max;
    }

    double I_mean_prev = emission.stats.Itot_mean;

    start_msg(t0, "=== Converging PIC with time window " + d2s(time_window, 2) + " fs\n");
    for (int i = 0; i < i_max; ++i) {
        solve_pic(time_window, i==0);
        emission.calc_global_stats();
        double err = (emission.stats.Itot_mean - I_mean_prev) / emission.stats.Itot_mean;
        if (MODES.VERBOSE){
            printf("  i=%d, I_mean= %e A, I_std=%.2f, error=%.2f\n", i, emission.stats.Itot_mean,
                    100. * emission.stats.Itot_std / emission.stats.Itot_mean, 100 * err);
        }
        I_mean_prev = emission.stats.Itot_mean;

        if (fabs(err) < 0.05 && fabs(err) < conf.SC.convergence * emission.stats.Itot_std /
                emission.stats.Itot_mean)
            return 0;
    }
    return 0;
}

int ProjectHeat::converge_heat(double T_ambient) {
    const int max_steps = 1000;
    double delta_time = conf.heating.delta_time;
    int ccg, hcg, step, error;

    bool global_verbosity = MODES.VERBOSE;

    start_msg(t0, "=== Converging heat...\n");

    for (step = 0; step < max_steps; ++step) {

//        emission.write("emission_before.movie");

        // advance heat and current system for delta_time
        error = solve_heat(conf.heating.t_ambient, delta_time, step == 0, ccg, hcg);
        if (error) return error;

        // modify the advance time depending on how slowly the solution is changing
        if (conf.field.mode == "laplace" || conf.field.mode == "converge")
            GLOBALS.TIME += delta_time;

        if (hcg < (ccg - 10) && delta_time <= conf.heating.dt_max / 1.25) // heat changed too little?
            delta_time *= 1.25;
        else if (hcg > (ccg + 10)) // heat changed too much?
            delta_time /= 1.25;

        // write debug data
        if (global_verbosity)
            printf( "t= %e ps, dt= %.2e ps, Tmax= %e K\n",
                    GLOBALS.TIME * 1.e-3, delta_time * 1.e-3, ch_solver.heat.max_solution() );
        write_results(true);

        // check if the result has converged
        if (max(hcg, ccg) < 10) break;

        // update field - advance PIC for delta time
        if (conf.field.mode == "transient")
            error = solve_pic(delta_time, false);
        else if (conf.field.mode == "converge")
            error = converge_pic(delta_time);
        if (error) return error;

    }

    MODES.VERBOSE = global_verbosity;
    end_msg(t0);

    check_return(step < max_steps, "Failed to converge heat equation after " + d2s(max_steps) + " steps!");
    return 0;
}

int ProjectHeat::write_results(bool force_write){

    if (!write_time() && !force_write) return 1;

    vacuum_interpolator.extract_solution(poisson_solver, conf.run.field_smoother);
    vacuum_interpolator.nodes.write("out/result_E_phi.movie");
    vacuum_interpolator.linhex.write("out/result_E_phi.vtk");

    if (conf.field.mode != "laplace"){
        emission.write("out/emission.movie");
        pic_solver.write("out/electrons.movie");
        surface_fields.write("out/surface_fields.movie");
        vacuum_interpolator.extract_charge_density(poisson_solver);
        vacuum_interpolator.nodes.write("out/result_E_charge.movie");
    }

    if (emission.atoms.size() > 0)
        emission.write("out/surface_emission.movie");

    if (conf.heating.mode != "none"){
        bulk_interpolator.nodes.write("out/result_J_T.movie");
        bulk_interpolator.lintet.write("out/result_J_T.vtk");
    }

    last_write_time = GLOBALS.TIME;
    return 0;
}

} /* namespace femocs */
