/*
 * ProjectHeat.cpp
 *
 *  Created on: 4.5.2018
 *      Author: kyritsak
 */

#include "ProjectSpaceCharge.h"

namespace femocs {

ProjectSpaceCharge::ProjectSpaceCharge(AtomReader &reader, Config &conf) : ProjectRunaway(reader, conf) {
}

int ProjectSpaceCharge::run(int timestep, double time) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (generate_mesh())
        return process_failed("Mesh generation failed!");

    check_return(!mesh_changed, "First meshing failed! Terminating...");

    cout << "Preparing solvers..." << endl;
    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");


    vector<double> I_target = {0.001};
    cout << find_Veff(I_target) << endl;
    //***** Run FEM solvers *****


//    for(auto factor : conf.field.apply_factors){
//        conf.field.E0 *= factor;
//        conf.field.V0 *= factor;
//
//        if (converge_pic())
//            return process_failed("Running field solver in a " + conf.field.solver + " mode failed!");
//
//        finalize(tstart);
//        conf.field.E0 /= factor;
//        conf.field.V0 /= factor;
//    }

    return 0;
}

int ProjectSpaceCharge::converge_pic() {
    double time_window = 32 * conf.pic.dt_max; //time window to check convergence
    int i_max = 1024; //window iterations

    double I_mean_prev = emission.global_data.I_mean;

    start_msg(t0, "=== Converging PIC...\n");
    for (int i = 0; i < i_max; ++i) {
        pic_solver.stats_reinit();
        solve_pic(time_window, i==0);
        emission.calc_global_stats();
        double err = (emission.global_data.I_mean - I_mean_prev) / emission.global_data.I_mean;
        if (MODES.VERBOSE){
            printf("  i=%d, I_mean= %e A, I_std=%.2f%, error=%.2f%, inj=%d, del=%d", i, emission.global_data.I_mean,
                    100. * emission.global_data.I_std / emission.global_data.I_mean, 100 * err, pic_solver.get_injected(), pic_solver.get_removed());
            cout << endl;
        }
        I_mean_prev = emission.global_data.I_mean;

        bool converged = fabs(err) < conf.pic.convergence * emission.global_data.I_std /
                emission.global_data.I_mean && fabs(err) < 0.05 && pic_solver.is_stable() ;

        if (converged)
            return 0;
    }
    return 0;
}

void ProjectSpaceCharge::get_currents(double Vappl, vector<double> &curs){
    curs.resize(conf.field.apply_factors.size());

    conf.emission.Vappl_SC = Vappl;

    int i = 0;
    for(auto factor : conf.field.apply_factors){
        solve_laplace(conf.field.E0 * factor, conf.field.V0 * factor);
        start_msg(t0, "=== Calculating electron emission...");

        cout << "interpolating surface fields" << endl;
        surface_fields.interpolate(ch_solver);

        cout << "inteprolating surface temperatures" << endl;
        if (!i) surface_temperatures.interpolate(ch_solver);

        cout << "initializing emission" << endl;
        emission.initialize(mesh, i == 0);

        cout << "calculating emission" << endl;
        emission.calc_emission(conf.emission, Vappl);
        end_msg(t0);
        curs[i++] = emission.global_data.I_tot;
    }

}

double ProjectSpaceCharge::find_Veff(vector<double> I_target){

    int Nmax = 20;
    double errlim = 0.01;

    bool lim_found = false;
    double Veff = conf.field.V0, Vhigh, Vlow;
    double errorlim;

    vector<double> currents;

    for(int i = 0; i < Nmax; ++i){
        get_currents(Veff, currents);
        double error = get_current_error(currents, I_target);

        cout << "Currents = " <<  currents[0] << " Veff = " << Veff << ", error = " <<  error << endl;
        if(error > errlim){
            Vhigh = Veff;
            if (lim_found)
                Veff = 0.5 * (Vhigh + Vlow);
            else
                Veff *= 2;
        } else if(error < -errlim) {
            Vlow = Veff;
            if (lim_found)
                Veff = 0.5 * (Vhigh + Vlow);
            else
                Veff /= 2;
        }
        else
            return Veff;
    }
}

double ProjectSpaceCharge::get_current_error(vector<double> I_calc, vector<double> I_target){
    require(I_calc.size() == I_target.size(), "comparison of current vectors no equal sizes");
    double error = 0;
    for(int i = 0; i < I_calc.size(); i++){
        error += log(I_calc[i] / I_target[i]);
    }
    return error;
}




} /* namespace femocs */
