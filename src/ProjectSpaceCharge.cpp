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
            find_Wsp();

            GLOBALS.TIME = 0;
            last_write_time = 0.;

            if (converge_pic())
                return process_failed("Running field solver in a " + conf.field.solver + " mode failed!");

            I_pic.push_back(emission.global_data.I_mean);

            export_on_line();
            string line_name = "out/line_data" + to_string((int) factor * 10) + ".dat";
            write_line(line_name);
        }
        conf.field.V0 = V_orig;
        conf.field.E0 = E_orig;
    }

    I_sc.resize(I_pic.size());

    double Veff = find_Veff();

    full_curve(Veff);



    write_output(Veff);

    return 0;
}

int ProjectSpaceCharge::converge_pic() {
    double time_window = 32 * conf.pic.dt_max; //time window to check convergence
//    conf.behaviour.write_period = time_window * (1. - numeric_limits<double>::epsilon());

    int i_max = 1024; //window iterations

    double I_mean_prev = emission.global_data.I_mean;

    start_msg(t0, "=== Converging PIC...\n");
    for (int i = 0; i < i_max; ++i) {
        pic_solver.stats_reinit();
        solve_pic(time_window, i == 0, true);
        emission.calc_global_stats();
        double err = (emission.global_data.I_mean - I_mean_prev) / emission.global_data.I_mean;

        if (MODES.VERBOSE){
            printf("  i=%d, I_mean= %e A, I_std=%.2f%, error=%.2f%, inj=%d, del=%d",
                    i, emission.global_data.I_mean,
                    100. * emission.global_data.I_std / emission.global_data.I_mean,
                    100 * err, pic_solver.get_injected(), pic_solver.get_removed());
            cout << endl;
        }

        I_mean_prev = emission.global_data.I_mean;

        bool converged = fabs(err) < conf.pic.convergence * emission.global_data.I_std /
                emission.global_data.I_mean && fabs(err) < 0.05;
        //&& pic_solver.is_stable() ;

        if (converged)
            return 0;
    }
    end_msg(t0);
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

void ProjectSpaceCharge::write_output(double Veff){
    ofstream out;
    out.open("out/results_SC.dat");
    out.setf(std::ios::scientific);
    out.precision(6);

    double Emax = surface_fields.calc_max_field();

    out << "effective Voltage = " << Veff << endl;
    out << "   F_max_L      Voltage       I_sc        I_pic" << endl;

    for (int i = 0; i < I_sc.size(); ++i)
        out << conf.field.apply_factors[i] * Emax << " " <<
                conf.field.apply_factors[i] * conf.field.V0 << " "
                << I_sc[i] << " " << I_pic[i] << endl;

    out.close();

    out.open("out/full_curve.dat");

    out << "   F_max_L      Voltage       I_sc        " << endl;

    for (int i = 0; i < I_full.size(); ++i)
        out << fact_full[i] * Emax << " " << fact_full[i] * conf.field.V0
                << " " << I_full[i] << endl;

    out.close();

}

void ProjectSpaceCharge::full_curve(double Veff){

    int Npoints = 128;
    double fmax = *max_element(conf.field.apply_factors.begin(), conf.field.apply_factors.end());
    double fmin = *min_element(conf.field.apply_factors.begin(), conf.field.apply_factors.end());

    fact_full.resize(Npoints);
    I_full.resize(Npoints);

    for(int i = 0; i < Npoints; ++i){
        fact_full[i] = fmin + i * (fmax - fmin) / (Npoints - 1);
        emission.set_sfactor(fact_full[i]);
        emission.calc_emission(conf.emission, Veff);
        I_full[i] = emission.global_data.I_tot;
    }
}

void ProjectSpaceCharge::find_Wsp(){
    double time_window = 16 * conf.pic.dt_max; //time window to check convergence
    int i_max = 1024; //window iterations

    double I_mean_prev = emission.global_data.I_mean;

    start_msg(t0, "=== Calculating a reasonable  Wsp ...\n");

    solve_pic(time_window, true);
    emission.calc_global_stats();
    double inj_per_step = pic_solver.get_injected() / (double) 32;
    if ((inj_per_step < 200 || inj_per_step > 1000))
        conf.pic.Wsp_el *= inj_per_step /  500;

    if (!inj_per_step){
        conf.pic.Wsp_el /= 1000;
        find_Wsp();
    }

    cout << "Found Wsp = " << conf.pic.Wsp_el << endl;
    pic_solver.reinit();

    end_msg(t0);
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
    for (int i = 0; i < N_points; i++)
        V_line[i] = fr.get_potential(i);

    vacuum_interpolator.extract_charge_density(poisson_solver);
    fr.calc_interpolation();

    for (int i = 0; i < N_points; i++)
        rho_line[i] = fr.get_potential(i);
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
