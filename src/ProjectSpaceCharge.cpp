/*
 * ProjectHeat.cpp
 *
 *  Created on: 4.5.2018
 *      Author: kyritsak
 */

#include "ProjectSpaceCharge.h"
#include <limits>

namespace femocs {

ProjectSpaceCharge::ProjectSpaceCharge(AtomReader &reader, Config &conf) : ProjectRunaway(reader, conf) {
    Ebase = conf.field.E0;
    Vbase = conf.field.V0;
    Fmax_base = Ebase;
}

int ProjectSpaceCharge::run(int timestep, double time) {
    double tstart = omp_get_wtime();

    if (generate_mesh())
        return process_failed("Mesh generation failed!");
    check_return(!mesh_changed, "First meshing failed! Terminating...");

    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    prepare_emission();

    if (conf.SC.I_pic.size()){
        cout << "Reading I_pic from file. Ipic = " << endl;
        for(int i = 0; i < conf.SC.I_pic.size(); ++i){
            I_pic.push_back(conf.SC.I_pic[i]);
            cout << I_pic[i] << endl;
        }

    } else{

        for(int i = 0; i < conf.SC.apply_factors.size(); ++i){
            double factor = conf.SC.apply_factors[i];
            prepare(i);
            GLOBALS.TIME = 0;

            if (converge_pic())
                return process_failed("Running field solver in a " + conf.field.mode + " mode failed!");

            /* TODO
             * What about replacing write_emission_stats with
             *   write_verbose_msg("fact=" + d2s(factor) + ", " + d2s(emission.stats));
             * */
            write_emission_stats("out/emission_stats_pic.dat", i == 0, factor);

            I_pic.push_back(emission.stats.Itot_mean);
//            export_on_line();
//            string line_name = "out/line_data" + to_string((int) (factor * 10.0)) + ".dat";
//            write_line(line_name);
        }
    }

    conf.field.V0 = Vbase;
    conf.field.E0 = Ebase;

    I_sc.resize(I_pic.size());

    double Veff = find_Veff();

    full_emission_curve(Veff);
    write_output(Veff);

    return 0;
}

int ProjectSpaceCharge::prepare(int i){

    int max_electrons = 500000;
    int max_Wsp_iter = 10;

    double init_factor, target_factor;

    start_msg(t0, "=== Preparing next applied field... \n");

    target_factor = conf.SC.apply_factors[i];
    if (i == 0)
        init_factor = target_factor / 10;
    else
        init_factor = conf.SC.apply_factors[i-1];

    if ((i == 0) || pic_solver.get_n_electrons() > max_electrons){
        pic_solver.reinit();
        init_factor = target_factor / 10;
    }

    for (int i = 0; i < max_Wsp_iter; i++){
        double inj_per_step = ramp_field(init_factor, target_factor);
        if (inj_per_step > 1.e100){
            pic_solver.reinit();
            conf.pic.weight_el *= 10;
            write_verbose_msg("Trying with higher Wsp");
        } else if (inj_per_step < 1){
            pic_solver.reinit();
            conf.pic.weight_el /= 10;
            write_verbose_msg("Trying with lower Wsp");
        } else if (inj_per_step < 200 || inj_per_step > 2000){
            pic_solver.reinit();
            conf.pic.weight_el *= inj_per_step /  600;
            init_factor = target_factor / 10; //reinitializing pic
            write_verbose_msg("Resetting the particle weight to Wsp = "
                    + to_string(conf.pic.weight_el));
        } else
            break;
    }
    end_msg(t0);
    return 0;
}

double ProjectSpaceCharge::ramp_field(double start_factor, double target_factor) {

    int window_steps = 12;
    double time_window = window_steps * conf.pic.dt_max; //time window to check convergence

    int ramping_steps = 5; // number of ramping window steps

    write_verbose_msg("ramping up the field from factor = " + to_string(start_factor) +
            "to " + to_string(target_factor));

    double factor_step = (target_factor - start_factor) / (ramping_steps - 1);

    for (int i = 0; i < ramping_steps; ++i) {
        double factor = start_factor + i * factor_step;
        conf.field.E0 = Ebase * factor;
        conf.field.V0 = Vbase * factor;
        pic_solver.stats_reinit(); //reinit stats to get the last cycle injections
        if (solve_pic(time_window, true))
            return 1.e200;
    }

    return pic_solver.get_injected() / (double) window_steps;
}

void ProjectSpaceCharge::prepare_emission() {
    solve_laplace(conf.field.E0, conf.field.V0);
    surface_fields.interpolate(ch_solver);
    surface_temperatures.interpolate(ch_solver);
    emission.initialize(mesh, true);

    emission.calc_emission(conf.emission, -1, true);

    Fmax_base = surface_fields.max_field();
}


int ProjectSpaceCharge::converge_pic() {

    int window_steps = 25;

    double time_window = window_steps * conf.pic.dt_max; //time window to check convergence

    int i_max = 32; //maximum window iterations to achieve convergence

    double I_mean_prev = emission.stats.Itot_mean;

    write_verbose_msg("Converging PIC...");

    for (int i = 0; i < i_max; ++i) {
        pic_solver.stats_reinit();
        solve_pic(time_window, i == 0);
        emission.calc_global_stats();
        double err = (emission.stats.Itot_mean - I_mean_prev) / emission.stats.Itot_mean;

        if (MODES.VERBOSE){
            printf("  i=%d, I_mean= %e A, I_std=%.2f, error=%.2f, inj=%d, del=%d",
                    i, emission.stats.Itot_mean,
                    100. * emission.stats.Itot_std / emission.stats.Itot_mean,
                    100 * err, pic_solver.get_injected(), pic_solver.get_removed());
            cout << endl;
        }

        I_mean_prev = emission.stats.Itot_mean;

        bool converged = fabs(err) < 5.e-3 || (fabs(err) < conf.SC.convergence * emission.stats.Itot_std /
                emission.stats.Itot_mean && fabs(err) < 0.05);
        //&& pic_solver.is_stable() ;

        if (converged)
            return 0;
    }
    return 0;
}

void ProjectSpaceCharge::get_currents(double Vappl){

    double Veff;
    for(int i = 0; i < I_pic.size(); ++i){
        emission.set_sfactor(conf.SC.apply_factors[i]);
        if (conf.emission.omega_SC > 0)
            Veff = Vappl * conf.SC.apply_factors[i];
        else
            Veff = Vappl;
        emission.calc_emission(conf.emission, Veff);
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

void ProjectSpaceCharge::write_output(double Veff){
    ofstream out;
    out.open("out/results_SC.dat");
    out.setf(std::ios::scientific);
    out.precision(6);

    if (conf.emission.omega_SC > 0)
        out << "omega_SC = " << Veff / conf.field.V0 << endl;
    else
        out << "effective Voltage = " << Veff << endl;

    out << "   F_max_L      Voltage       I_sc        I_pic" << endl;

    for (int i = 0; i < I_sc.size(); ++i)
        out << conf.SC.apply_factors[i] * Fmax_base << " " <<
                conf.SC.apply_factors[i] * conf.field.V0 << " "
                << I_sc[i] << " " << I_pic[i] << endl;

    out.close();

}

void ProjectSpaceCharge::full_emission_curve(double Vappl){

    int Npoints = 128;
    double fmax = *max_element(conf.SC.apply_factors.begin(), conf.SC.apply_factors.end());
    double fmin = *min_element(conf.SC.apply_factors.begin(), conf.SC.apply_factors.end());
    double factor, Veff;

    for(int i = 0; i < Npoints; ++i){
        factor = fmin + i * (fmax - fmin) / (Npoints - 1);

        if (conf.emission.omega_SC > 0)
            Veff = Vappl * factor;
        else
            Veff = Vappl;

        emission.set_sfactor(factor);
        emission.calc_emission(conf.emission, Veff);
        write_emission_data("out/emission_full_data.dat", i == 0);
    }
}

void ProjectSpaceCharge::export_on_line(){

    const int interpolation_rank = 3;
    const int N_points = 2048;

    z_line.resize(N_points);
    V_line.resize(N_points);
    rho_line.resize(N_points);

    FieldReader fr(&vacuum_interpolator);
    fr.set_preferences(false, 3, interpolation_rank);
    fr.reserve(N_points);

    if (!apex_found) find_apex();

    for (int i = 0; i < N_points; i++){
        z_line[i] = exp(log(apex.z) + i * (log(anode.z) - log(apex.z)) / (N_points - 1));
        Point3 point = Point3(apex.x, apex.y, z_line[i]);
        fr.append(point);
    }

    vacuum_interpolator.extract_solution(poisson_solver, false);
    fr.calc_interpolation();

    fr.calc_interpolation();
    for (int i = 0; i < N_points; i++) {
        V_line[i] = fr.get_potential(i);
        rho_line[i] = fr.get_charge_dens(i);
    }
}


void ProjectSpaceCharge::write_line(string filename){
    ofstream out;
    out.open(filename);
    out.setf(std::ios::scientific);
    out.precision(6);

    out << "     z       Potential       Charge density  " << endl;

    for (int i = 0; i < V_line.size(); ++i)
        out << z_line[i] - apex.z << " " << V_line[i] << " " << rho_line[i] << endl;

    out.close();
}

void ProjectSpaceCharge::write_emission_stats(string filename, bool first_time, double factor){
    ofstream out;
    out.open(filename, ios_base::app);
    out.setf(std::ios::scientific);
    out.precision(6);

    if (first_time) {
        out << "factor      ";
        emission.write_stats(out, true);
    }

    out << factor << " ";
    emission.write_stats(out, false);
    out.close();
}

void ProjectSpaceCharge::write_emission_data(string filename, bool first_time){
    ofstream out;
    out.open(filename, ios_base::app);
    out.setf(std::ios::scientific);
    out.precision(6);

    if (first_time) {
        out << "Fmax_base = " << Fmax_base << "\n";
        emission.write_global_data(out, true);
    }

    emission.write_global_data(out, false);
    out.close();
}



void ProjectSpaceCharge::find_apex(){
    int i_max;
    double zmax = -1.e-100;

    for (int i = 0; i < surface_fields.size(); i++){
        if (surface_fields.get_point(i).z > zmax){
            i_max = i;
            zmax = surface_fields.get_point(i).z;
        }
    }
    apex = Point3(surface_fields.get_point(i_max));

    for (int i = 0; i < vacuum_interpolator.nodes.size(); i++){
        if (vacuum_interpolator.nodes.get_vertex(i).z > zmax){
            i_max = i;
            zmax = vacuum_interpolator.nodes.get_vertex(i).z;
        }
    }
    anode = Point3(vacuum_interpolator.nodes.get_vertex(i_max));
    apex_found = true;
}


} /* namespace femocs */
