/*
 * ProjectHeat.cpp
 *
 *  Created on: 4.5.2018
 *      Author: kyritsak
 */

#include "ProjectHeat.h"

namespace femocs {

ProjectHeat::ProjectHeat(AtomReader &reader, Config &conf) : ProjectSpaceCharge(reader, conf){
    last_heat_time = 0.0;
}

int ProjectHeat::run(int timestep, double time) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (reinit()) {
        write_verbose_msg("Atoms haven't moved significantly. Previous mesh will be used.");
        dense_surf.update_positions(reader);
    }
    else if (generate_mesh())
        return process_failed("Mesh generation failed!");

    check_return(!mesh_changed, "First meshing failed! Terminating...");

    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    //***** Run FEM solvers *****

    for(auto factor : conf.scharge.apply_factors){
        conf.field.E0 *= factor;
        conf.field.V0 *= factor;

        if (run_field_solver())
            return process_failed("Running field solver in a " + conf.field.mode + " mode failed!");

        if (run_heat_solver())
            return process_failed("Running heat solver in a " + conf.heating.mode + " mode failed!");

        update_mesh_pointers();

        //***** Prepare for data export and next run *****
        if (prepare_export())
            return process_failed("Interpolating solution on atoms failed!");

        finalize(tstart);
        conf.field.E0 /= factor;
        conf.field.V0 /= factor;
    }

    return 0;
}

int ProjectHeat::run_field_solver() {
    if (conf.field.mode == "transient")
        return solve_pic(conf.behaviour.timestep_fs, mesh_changed);
    else if (conf.field.mode == "converge")
        return converge_pic();
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

int ProjectHeat::converge_heat(double T_ambient) {
    const int max_steps = 1000;
    double delta_time = conf.heating.delta_time;
    int ccg, hcg, step, error;

    for (step = 0; step < max_steps; ++step) {

        GLOBALS.TIME = last_heat_time + delta_time;
        // advance heat and current system for delta_time
        error = solve_heat(conf.heating.t_ambient, delta_time, step == 0, ccg, hcg);
        if (error) return error;

        if (MODES.VERBOSE)
            printf( "t= %g ps, dt= %.2g ps \n",
                    GLOBALS.TIME * 1.e-3, delta_time * 1.e-3);

        vacuum_interpolator.nodes.write("out/result_E_phi.movie", true);
        bulk_interpolator.nodes.write("out/result_J_T.movie", true);

        if (hcg < (ccg - 10) && delta_time <= conf.heating.dt_max / 1.25) // heat changed too little?
            delta_time *= 1.25;
        else if (hcg > (ccg + 10)) // heat changed too much?
            delta_time /= 1.25;

        // check if the result has converged
        if (hcg < 10) break;

        calc_surf_temperatures();

        // update field - advance PIC for delta time
        if (conf.field.mode == "transient")
            error = solve_pic(delta_time, false);
        else if (conf.field.mode == "converge")
            error = converge_pic();
        if (error) return error;

    }

    check_return(step < max_steps, "Failed to converge heat equation after " + d2s(max_steps) + " steps!");
    return 0;
}

} /* namespace femocs */
