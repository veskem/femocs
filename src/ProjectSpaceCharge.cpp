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

    if (generate_mesh())
        return process_failed("Mesh generation failed!");
    check_return(!mesh_changed, "First meshing failed! Terminating...");

    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");


    if (conf.emission.I_pic.size()){
        for(int i = 0; i < conf.emission.I_pic.size(); ++i)
            I_pic.push_back(conf.emission.I_pic[i]);
    } else{
        double E_orig = conf.field.E0, V_orig = conf.field.V0;
        for(auto factor : conf.field.apply_factors){
            conf.field.E0 = E_orig * factor;
            conf.field.V0 = V_orig * factor;

            if (converge_pic())
                return process_failed("Running field solver in a " + conf.field.solver + " mode failed!");

            I_pic.push_back(emission.global_data.I_mean);
        }
        conf.field.V0 = V_orig;
        conf.field.E0 = E_orig;
    }

    I_sc.resize(I_pic.size());

    double Veff = find_Veff();

    write_results(Veff);

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

void ProjectSpaceCharge::get_currents(double Vappl){

    for(int i = 0; i < I_pic.size(); ++i){
        emission.set_sfactor(conf.field.apply_factors[i]);
        emission.calc_emission(conf.emission, Vappl);
        I_sc[i] = emission.global_data.I_tot;
    }
}

double ProjectSpaceCharge::find_Veff(){


    solve_laplace(conf.field.E0, conf.field.V0);
    surface_fields.interpolate(ch_solver);
    surface_temperatures.interpolate(ch_solver);
    emission.initialize(mesh, true);

    int Nmax = 50;
    double errlim = 0.01;

    double Veff = conf.field.V0, Vhigh, Vlow, old_error;

    vector<double> currents;

    for(int i = 0; i < Nmax; ++i){ // first find two values that produce opposite sign errors
        get_currents(Veff);
        double error = get_current_error();
        if (i == 0)
            old_error = error;
        cout << " Veff = " << Veff << ", error = " <<  error << endl;
        if(error > errlim && error * old_error > 0){
            Vlow = Veff;
            Veff *= 2;
            Vhigh = Veff;
        } else if(error < -errlim && error * old_error > 0){
            Vhigh = Veff;
            Veff /= 2;
            Vlow = Veff;
        } else
            break;
        old_error = error;
    }

    Veff = .5 * (Vhigh + Vlow);

    for(int i = 0; i < Nmax; ++i){ // perform bisection
        get_currents(Veff);
        double error = get_current_error();
        cout  << " Veff = " << Veff << ", error = " <<  error << endl;
        if(error > errlim){
            Vlow = Veff;
            Veff = .5 * (Vhigh + Vlow);
        }
        else if(error < -errlim){
            Vhigh = Veff;
            Veff = .5 * (Vhigh + Vlow);
        }
        else
            break;
    }

    return Veff;
}

double ProjectSpaceCharge::get_current_error(){
    require(I_sc.size() == I_pic.size(), "comparison of current vectors no equal sizes");
    double error = 0;
    for(int i = 0; i < I_sc.size(); i++){
        error += log(I_sc[i] / I_pic[i]);
    }
    return error;
}

void ProjectSpaceCharge::write_results(double Veff){
    ofstream out;
    out.open("results_SC.dat");
    out.setf(std::ios::scientific);
    out.precision(6);

    double Emax = surface_fields.calc_max_field();

    out << "effective Voltage = " << Veff << endl;
    out << "   F_max_L      Voltage       I_sc        I_pic" << endl;

    for (int i = 0; i < I_sc.size(); ++i){
        out << conf.field.apply_factors[i] * Emax << " " <<
                conf.field.apply_factors[i] * conf.field.V0 << " "
                << I_sc[i] << " " << I_pic[i] << endl;
    }

    out.close();

}


} /* namespace femocs */
