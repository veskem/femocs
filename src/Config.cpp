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

// Config constructor initializes configuration parameters
Config::Config() {
    path.extended_atoms = "";         // file with the atoms forming the extended surface
    path.infile = "";                 // file with the nanostructure atoms
    path.mesh_file = "";              // file containing triangular and tetrahedral mesh data

    behaviour.verbosity = "verbose";  // mute, silent, verbose
    behaviour.project = "runaway";    // project to be run; runaway, ...
    behaviour.n_writefile = 1;        // number of time steps between writing the output files
    behaviour.interpolation_rank = 1; // rank of the solution interpolation; 1-linear tetrahedral, 2-quadratic tetrahedral, 3-linear hexahedral
    behaviour.write_period = 1.e5;    // write files every write_period of time
    behaviour.timestep_fs = 4.05;     // Total time of a FEMOCS run [fs]
    behaviour.mass = 63.5460;         // Atom mass [amu]
    behaviour.rnd_seed = 12345;       // Seed for random number generator
    behaviour.n_omp_threads = 1;      // Number of opened OpenMP threads

    run.cluster_anal = true;          // enable cluster analysis
    run.apex_refiner = false;         // refine nanotip apex
    run.rdf = false;                  // use radial distribution function to recalculate lattice constant and coordination analysis parameters
    run.output_cleaner = true;        // clear output folder
    run.surface_cleaner = true;       // clean surface by measuring the atom distance from the triangular surface
    run.field_smoother = true;        // replace nodal field with the average of its neighbouring nodal fields

    geometry.mesh_quality = "2.0";    // minimum tetrahedron quality Tetgen is allowed to make
    geometry.element_volume = "";     // maximum tetrahedron volume Tetgen is allowed to make
    geometry.nnn = 12;                // number of nearest neighbours in bulk
    geometry.latconst = 3.61;         // lattice constant
    geometry.coordination_cutoff = 3.1; // coordination analysis cut-off radius
    geometry.cluster_cutoff = 0;      // cluster analysis cut-off radius
    geometry.charge_cutoff = 30;      // Coulomb force cut-off radius
    geometry.surface_thickness = 3.1; // maximum distance the surface atom is allowed to be from surface mesh
    geometry.box_width = 10;          // minimal simulation box width in units of tip height
    geometry.box_height = 6;          // simulation box height in units of tip height
    geometry.bulk_height = 20;        // bulk substrate height [lattice constant]
    geometry.radius = 0.0;            // inner radius of coarsening cylinder
    geometry.height = 0.0;            // height of generated artificial nanotip in the units of radius

    tolerance.charge_min = 0.8;       // min ratio face charges are allowed to deviate from the total charge
    tolerance.charge_max = 1.2;       // max ratio face charges are allowed to deviate from the total charge
    tolerance.field_min = 0.1;        // min ratio numerical field can deviate from analytical one
    tolerance.field_max = 5.0;        // max ratio numerical field can deviate from analytical one
    tolerance.distance = 0.0;         // rms distance tolerance for atom movement between two time steps

    field.E0 = 0.0;                   // long range electric field
    field.ssor_param = 1.2;           // parameter for SSOR preconditioner
    field.phi_error = 1e-9;           // maximum allowed electric potential error
    field.n_phi = 10000;              // maximum number of Conjugate Gradient iterations in phi calculation
    field.V0 = 0.0;                   // anode voltage
    field.anode_BC = "neumann";        // anode Neumann boundary
    field.solver = "laplace";         // type of field equation to be solved; laplace or poisson
    field.element_degree = 1;         // FEM element shape function degree

    heating.mode = "none";            // method to calculate current density and temperature; none, stationary or transient
    heating.rhofile = "in/rho_table.dat"; // rho table file
    heating.lorentz = 2.44e-8;        // Lorentz number
    heating.t_ambient = 300.0;        // ambient temperature
    heating.t_error = 10.0;           // max allowed temperature error in Newton iterations
    heating.n_newton = 10;            // max number of Newton iterations
    heating.n_cg = 2000;              // max number of Conjugate-Gradient iterations
    heating.cg_tolerance = 1e-9;      // solution accuracy in Conjugate-Gradient solver
    heating.ssor_param = 1.2;         // parameter for SSOR pre-conditioner in DealII; 1.2 is known to work well with Laplace
    heating.delta_time = 1.e3;        // timestep of time domain integration [fs]
    heating.dt_max = 1.e5;            // max allowed timestep
    heating.tau = 100.0;              // time constant in Berendsen thermostat
    heating.assemble_method = "euler"; // method to assemble system matrix for solving heat equation

    emission.blunt = true;            // by default emitter is blunt (simple SN barrier used for emission)
    emission.cold = false;
    emission.work_function = 4.5;     // work function [eV]
    emission.omega_SC = -1;           // SC is ignored in Emission by default
    emission.SC_error = 1.e-3;        // Convergence criterion for SC iteration
    emission.Vappl_SC = 0.;           // Vappl used for SC calculations
    emission.SC_mode = "local";       // local current density with SC

    force.mode = "none";              // forces to be calculated; lorentz, all, none

    smoothing.algorithm = "laplace";  // surface mesh smoother algorithm; none, laplace or fujiwara
    smoothing.n_steps = 0;            // number of surface mesh smoothing iterations
    smoothing.lambda_mesh = 0.6307;   // lambda parameter in surface mesh smoother
    smoothing.mu_mesh = -0.6732;      // mu parameter in surface mesh smoother
    smoothing.beta_atoms = 0.0;       // surface smoothing factor; bigger number gives smoother surface
    smoothing.beta_charge = 1.0;      // charge smoothing factor; bigger number gives smoother charges

    cfactor.amplitude = 0.4;          // coarsening factor outside the warm region
    cfactor.r0_cylinder = 0;          // minimum distance between atoms in nanotip outside the apex
    cfactor.r0_sphere = 0;            // minimum distance between atoms in nanotip apex
    cfactor.exponential = 0.5;        // coarsening rate; min distance between coarsened atoms outside the warm region is
                                      // d_min ~ pow(|r1-r2|, exponential)
    pic.mode = "none";
    pic.dt_max = 1.0;
    pic.Wsp_el =  .01;
    pic.fractional_push = true;
    pic.coll_coulomb_ee = false;

    SC.convergence = .1;
}

// Remove the noise from the beginning of the string
void Config::trim(string& str) {
    str.erase(0, str.find_first_of(comment_symbols + data_symbols));
}

// Read the configuration parameters from input script
void Config::read_all(const string& file_name) {
    if(file_name == "") return;

    // Store the commands and their arguments
    parse_file(file_name);

    // Check for the obsolete commands
    check_obsolete("postprocess_marking");
    check_obsolete("force_factor");
    check_obsolete("heating", "heating_mode");
    check_obsolete("surface_thichness", "surface_thickness");
    check_obsolete("smooth_factor", "surface_smooth_factor");
    check_obsolete("surface_cleaner", "clean_surface");
    check_obsolete("run_pic", "pic_mode");
    check_obsolete("use_histclean");

    // Modify the parameters that are specified in input script
    read_command("work_function", emission.work_function);
    read_command("emitter_blunt", emission.blunt);
    read_command("omega_SC", emission.omega_SC);
    read_command("maxerr_SC", emission.SC_error);
    read_command("emitter_cold", emission.cold);
    read_command("Vappl_SC", emission.Vappl_SC);
    read_command("SC_mode", emission.SC_mode);

    read_command("t_ambient", heating.t_ambient);
    read_command("heating_mode", heating.mode);
    read_command("lorentz", heating.lorentz);
    read_command("rhofile", heating.rhofile);
    read_command("heat_dtinit", heating.delta_time);
    read_command("heat_dtmax", heating.dt_max);

    read_command("smooth_steps", smoothing.n_steps);
    read_command("smooth_lambda", smoothing.lambda_mesh);
    read_command("smooth_mu", smoothing.mu_mesh);
    read_command("smooth_algorithm", smoothing.algorithm);
    read_command("surface_smooth_factor", smoothing.beta_atoms);
    read_command("charge_smooth_factor", smoothing.beta_charge);

    read_command("phi_error", field.phi_error);
    read_command("n_phi", field.n_phi);
    read_command("elfield", field.E0);
    read_command("Vappl", field.V0);
    read_command("anode_BC", field.anode_BC);
    read_command("field_solver", field.solver);
    read_command("element_degree", field.element_degree);

    read_command("force_mode", force.mode);

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
    read_command("write_log", MODES.WRITELOG);

    read_command("femocs_verbose_mode", behaviour.verbosity);
    read_command("project", behaviour.project);
    read_command("n_writefile", behaviour.n_writefile);
    read_command("interpolation_rank", behaviour.interpolation_rank);
    read_command("write_period", behaviour.write_period);
    read_command("femocs_run_time", behaviour.timestep_fs);
    read_command("mass(1)", behaviour.mass);
    read_command("seed", behaviour.rnd_seed);
    read_command("n_omp", behaviour.n_omp_threads);

    read_command("distance_tol", tolerance.distance);

    read_command("pic_mode", pic.mode);
    read_command("pic_dtmax", pic.dt_max);
    read_command("electronWsp", pic.Wsp_el);
    read_command("pic_fractional_push", pic.fractional_push);
    read_command("pic_collide_coulomb_ee", pic.coll_coulomb_ee);
    
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

// Read the commands and their arguments from the file and store them into the buffer
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

// Check for the obsolete commands from the buffered commands
void Config::check_obsolete(const string& command) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            write_verbose_msg("Command '" + command + "' is obsolete! You can safely remove it!");
            return;
        }
}

// Check for the obsolete commands that are similar to valid ones
void Config::check_obsolete(const string& command, const string& substitute) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            write_verbose_msg("Command '" + command + "' is obsolete!"
                    " It is similar yet different to the command '" + substitute + "'!");
            return;
        }
}

// Look up the parameter with string argument
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

// Look up the parameter with boolean argument
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

// Look up the parameter with unsigned integer argument
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

// Look up the parameter with integer argument
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

// Look up the parameter with double argument
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

// Print the stored commands and parameters
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
