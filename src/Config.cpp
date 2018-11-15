/*
 * Config.cpp
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#include <fstream>
#include <algorithm>

#include "Config.h"
#include "Globals.h"

using namespace std;
namespace femocs {

Config::Config() {
    path.extended_atoms = "";
    path.infile = "";
    path.mesh_file = "";

    behaviour.verbosity = "verbose";
    behaviour.project = "runaway";
    behaviour.n_write_log = -1;
    behaviour.n_writefile = 1;
    behaviour.interpolation_rank = 1;
    behaviour.n_read_conf = 0;
    behaviour.write_period = 4.05;
    behaviour.timestep_fs = 4.05;
    behaviour.mass = 63.5460;
    behaviour.rnd_seed = 12345;
    behaviour.n_omp_threads = 1;
    behaviour.timestep_step = 1;

    run.cluster_anal = true;
    run.apex_refiner = false;
    run.rdf = false;
    run.output_cleaner = true;
    run.surface_cleaner = true;
    run.field_smoother = true;

    geometry.mesh_quality = "2.0";
    geometry.element_volume = "";
    geometry.nnn = 12;
    geometry.latconst = 3.61;
    geometry.coordination_cutoff = 3.1;
    geometry.cluster_cutoff = 0;
    geometry.charge_cutoff = 30;
    geometry.surface_thickness = 3.1;
    geometry.box_width = 10;
    geometry.box_height = 6;
    geometry.bulk_height = 20;
    geometry.radius = 0.0;
    geometry.height = 0.0;
    geometry.distance_tol = 0.0;

    tolerance.charge_min = 0.8;
    tolerance.charge_max = 1.2;
    tolerance.field_min = 0.1;
    tolerance.field_max = 5.0;

    field.E0 = 0.0;
    field.ssor_param = 1.2;
    field.cg_tolerance = 1e-9;
    field.n_cg = 10000;
    field.V0 = 0.0;
    field.anode_BC = "neumann";
    field.mode = "laplace";
    field.assemble_method = "parallel";

    heating.mode = "none";
    heating.rhofile = "in/rho_table.dat";
    heating.lorentz = 2.44e-8;
    heating.t_ambient = 300.0;
    heating.n_cg = 2000;
    heating.cg_tolerance = 1e-9;
    heating.ssor_param = 1.2;         // 1.2 is known to work well with Laplace
    heating.delta_time = 10.0;
    heating.dt_max = 1.0e5;
    heating.tau = 100.0;
    heating.assemble_method = "euler";

    emission.blunt = true;
    emission.cold = false;
    emission.work_function = 4.5;
    emission.omega_SC = -1;
    emission.Vappl_SC = 0.;

    force.mode = "none";

    smoothing.algorithm = "laplace";
    smoothing.n_steps = 0;
    smoothing.lambda_mesh = 0.6307;
    smoothing.mu_mesh = -0.6732;
    smoothing.beta_atoms = 0.0;
    smoothing.beta_charge = 1.0;

    cfactor.amplitude = 0.4;
    cfactor.r0_cylinder = 0;
    cfactor.r0_sphere = 0;
    cfactor.exponential = 0.5;

    pic.dt_max = 1.0;
    pic.weight_el =  .01;
    pic.fractional_push = true;
    pic.collide_ee = true;
    pic.periodic = false;
    pic.landau_log = 13.0;

    SC.convergence = .1;
}

void Config::trim(string& str) {
    str.erase(0, str.find_first_of(comment_symbols + data_symbols));
}

void Config::read_all() {
    read_all(file_name, false);
}

void Config::read_all(const string& fname, bool full_run) {
    if (fname == "") return;
    file_name = fname;

    // Store the commands and their arguments
    parse_file(fname);

    if (full_run) {
        // Check for the obsolete commands
        check_obsolete("postprocess_marking");
        check_obsolete("force_factor");
        check_obsolete("use_histclean");
        check_obsolete("maxerr_SC");
        check_obsolete("SC_mode");

        // Check for the changed commands
        check_changed("heating", "heat_mode");
        check_changed("surface_thichness", "surface_thickness");
        check_changed("smooth_factor", "surface_smooth_factor");
        check_changed("surface_cleaner", "clean_surface");
        check_changed("write_log", "n_write_log");
        check_changed("run_pic", "pic_mode");
        check_changed("electronWsp", "electron_weight");
        check_changed("field_solver", "field_mode");
        check_changed("pic_mode", "field_mode");
        check_changed("heating_mode", "heat_mode");
    }

    // Modify the parameters that are specified in input script
    read_command("work_function", emission.work_function);
    read_command("emitter_blunt", emission.blunt);
    read_command("omega_SC", emission.omega_SC);
    read_command("emitter_cold", emission.cold);
    read_command("Vappl_SC", emission.Vappl_SC);

    read_command("heat_mode", heating.mode);
    read_command("rhofile", heating.rhofile);
    read_command("lorentz", heating.lorentz);
    read_command("t_ambient", heating.t_ambient);
    read_command("heat_ncg", heating.n_cg);
    read_command("heat_cgtol", heating.cg_tolerance);
    read_command("heat_ssor", heating.ssor_param);
    read_command("heat_dt", heating.delta_time);
    read_command("heat_dtmax", heating.dt_max);
    read_command("vscale_tau", heating.tau);
    read_command("heat_assemble", heating.assemble_method);

    read_command("field_mode", field.mode);
    read_command("field_ssor", field.ssor_param);
    read_command("field_cgtol", field.cg_tolerance);
    read_command("field_ncg", field.n_cg);
    read_command("elfield", field.E0);
    read_command("Vappl", field.V0);
    read_command("anode_BC", field.anode_BC);
    read_command("field_assemble", field.assemble_method);

    read_command("force_mode", force.mode);

    read_command("smooth_steps", smoothing.n_steps);
    read_command("smooth_lambda", smoothing.lambda_mesh);
    read_command("smooth_mu", smoothing.mu_mesh);
    read_command("smooth_algorithm", smoothing.algorithm);
    read_command("surface_smooth_factor", smoothing.beta_atoms);
    read_command("charge_smooth_factor", smoothing.beta_charge);

    read_command("force_mode", force.mode);

    read_command("distance_tol", geometry.distance_tol);
    read_command("latconst", geometry.latconst);
    read_command("coord_cutoff", geometry.coordination_cutoff);
    read_command("cluster_cutoff", geometry.cluster_cutoff);
    read_command("charge_cutoff", geometry.charge_cutoff);
    read_command("surface_thickness", geometry.surface_thickness);
    read_command("nnn", geometry.nnn);
    read_command("mesh_quality", geometry.mesh_quality);
    read_command("element_volume", geometry.element_volume);
    read_command("radius", geometry.radius);
    read_command("tip_height", geometry.height);
    read_command("box_width", geometry.box_width);
    read_command("box_height", geometry.box_height);
    read_command("bulk_height", geometry.bulk_height);

    read_command("extended_atoms", path.extended_atoms);
    read_command("infile", path.infile);
    read_command("mesh_file", path.mesh_file);

    read_command("cluster_anal", run.cluster_anal);
    read_command("refine_apex", run.apex_refiner);
    read_command("use_rdf", run.rdf);
    read_command("clear_output", run.output_cleaner);
    read_command("clean_surface", run.surface_cleaner);
    read_command("smoothen_field", run.field_smoother);
    read_command("femocs_periodic", MODES.PERIODIC);

    read_command("n_read_conf", behaviour.n_read_conf);
    read_command("n_write_log", behaviour.n_write_log);
    read_command("femocs_verbose_mode", behaviour.verbosity);
    read_command("project", behaviour.project);
    read_command("n_writefile", behaviour.n_writefile);
    read_command("interpolation_rank", behaviour.interpolation_rank);
    read_command("write_period", behaviour.write_period);
    read_command("md_timestep", behaviour.timestep_fs);
    read_command("mass(1)", behaviour.mass);
    read_command("seed", behaviour.rnd_seed);
    read_command("n_omp", behaviour.n_omp_threads);
    read_command("timestep_step", behaviour.timestep_step);

    read_command("pic_dtmax", pic.dt_max);
    read_command("electron_weight", pic.weight_el);
    read_command("pic_fractional_push", pic.fractional_push);
    read_command("pic_collide_ee", pic.collide_ee);
    read_command("pic_periodic", pic.periodic);
    read_command("pic_landau_log", pic.landau_log);
    
    read_command("SC_converge_criterion", SC.convergence);

    // Read commands with potentially multiple arguments like...
    vector<double> args;
    int n_read_args;

    // ...charge and field tolerances
    args = {0, 0};
    n_read_args = read_command("charge_tolerance", args);
    if (n_read_args == 1) {
        tolerance.charge_min = 1.0 - args[0];
        tolerance.charge_max = 1.0 + args[0];
    } else if (n_read_args == 2) {
        tolerance.charge_min = args[0];
        tolerance.charge_max = args[1];
    }

    n_read_args = read_command("field_tolerance", args);
    if (n_read_args == 1) {
        tolerance.field_min = 1.0 - args[0];
        tolerance.field_max = 1.0 + args[0];
    } else if (n_read_args == 2) {
        tolerance.field_min = args[0];
        tolerance.field_max = args[1];
    }

    // ...coarsening factors
    read_command("coarse_rate", cfactor.exponential);
    args = {cfactor.amplitude, (double)cfactor.r0_cylinder, (double)cfactor.r0_sphere};
    n_read_args = read_command("coarse_factor", args);
    cfactor.amplitude = args[0];
    cfactor.r0_cylinder = static_cast<int>(args[1]);
    cfactor.r0_sphere = static_cast<int>(args[2]);

    SC.apply_factors.resize(128);
    n_read_args = read_command("apply_factors", SC.apply_factors);
    if (n_read_args > 0)
        SC.apply_factors.resize(n_read_args);
    else
        SC.apply_factors = {1.};

    SC.I_pic.resize(128);
    n_read_args = read_command("currents_pic", SC.I_pic);
    SC.I_pic.resize(n_read_args);
    require(!n_read_args ||n_read_args == field.apply_factors.size(),
            "current_pic & apply_factors sizes don't match");

}

void Config::parse_file(const string& file_name) {
    ifstream file(file_name);
    require(file, "File not found: " + file_name);

    string line;
    data.clear();

    // loop through the lines in a file
    while (getline(file, line)) {
        line += " "; // needed to find the end of line

        bool line_started = true;
        // store the command and its parameters from non-empty and non-pure-comment lines
        while(line.size() > 0) {
            trim(line);
            int i = line.find_first_not_of(data_symbols);
            if (i <= 0) break;

            if (line_started && line.substr(0, i) == "femocs_end") return;
            if (line_started) data.push_back({});
            if (line_started) {
                // force all the characters in a command to lower case
                string command = line.substr(0, i);
                std::transform(command.begin(), command.end(), command.begin(), ::tolower);
                data.back().push_back(command);
            } else {
                data.back().push_back( line.substr(0, i) );
            }

            line = line.substr(i);
            line_started = false;
        }
    }
}

void Config::check_obsolete(const string& command) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            write_silent_msg("Command '" + command + "' is obsolete! You can safely remove it!");
            return;
        }
}

void Config::check_changed(const string& command, const string& substitute) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            write_silent_msg("Command '" + command + "' has changed!"
                    " It is similar yet different to the command '" + substitute + "'!");
            return;
        }
}

int Config::read_command(string param, string& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            arg = str[1]; return 0;
        }
    return 1;
}

int Config::read_command(string param, bool& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is1(str[1]);
            istringstream is2(str[1]);
            bool result;
            // try to parse the bool argument in text format
            if (is1 >> std::boolalpha >> result) { arg = result; return 0; }
            // try to parse the bool argument in numeric format
            else if (is2 >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, unsigned int& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); int result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, int& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); int result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, double& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); double result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, vector<double>& args) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    int n_read_args = 0;
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param)
            for (unsigned i = 0; i < args.size() && i < (str.size()-1); ++i) {
                istringstream is(str[i+1]);
                double result;
                if (is >> result) { args[i] = result; n_read_args++; }
            }
    return n_read_args;
}

void Config::print_data() {
    if (!MODES.VERBOSE) return;
    const int cmd_len = 20;

    for (const vector<string>& line : data) {
        for (const string& ln : line) {
            int str_len = ln.length();
            int whitespace_len = max(1, cmd_len - str_len);
            cout << ln << string(whitespace_len, ' ');
        }

        cout << endl;
    }
}

} // namespace femocs
