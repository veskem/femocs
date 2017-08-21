/*
 * Config.cpp
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#include "Config.h"
#include <fstream>
#include <algorithm>

using namespace std;
namespace femocs {

// Config constructor initializes configuration parameters
Config::Config() : obsolete_commands{
    "postprocess_marking",
    "heating",
    "force_factor"
} {

    extended_atoms = "";         // file with the atoms forming the extended surface
    atom_file = "";              // file with the nanostructure atoms
    latconst = 3.61;             // lattice constant
    surface_thichness = 3.1;     // maximum distance the surface atom is allowed to be from surface mesh
    coordination_cutoff = 3.1;   // coordination analysis cut-off radius
    cluster_cutoff = 0;          // cluster analysis cut-off radius
    nnn = 12;                    // number of nearest neighbours in bulk
    mesh_quality = "2.0";        // minimum tetrahedron quality Tetgen is allowed to make
    element_volume = "";         // maximum tetrahedron volume Tetgen is allowed to make
    radius = 0.0;                // inner radius of coarsening cylinder
    box_width = 10;              // minimal simulation box width in units of tip height
    box_height = 6;              // simulation box height in units of tip height
    bulk_height = 20;            // bulk substrate height [lattice constant]

    surface_smooth_factor = 0.0; // surface smoothing factor; bigger number gives smoother surface
    charge_smooth_factor = 1.0;  // charge smoothing factor; bigger number gives smoother charges
    cfactor.amplitude = 0.4;     // coarsening factor outside the warm region
    cfactor.r0_cylinder = 0;     // minimum distance between atoms in nanotip outside the apex
    cfactor.r0_sphere = 0;       // minimum distance between atoms in nanotip apex
    t_ambient = 300.0;           // ambient temperature
    t_error = 10.0;              // maximum allowed temperature error in Newton iterations
    n_newton = 10;               // maximum number of Newton iterations
    phi_error = 1e-9;            // maximum allowed electric potential error
    n_phi = 10000;               // maximum number of Conjugate Gradient iterations in phi calculation
    ssor_param = 1.2;            // parameter for SSOR preconditioner
    charge_tolerance = 0.2;      // tolerance how much face charges are allowed to deviate from the long range one
    refine_apex = false;         // refine nanotip apex
    distance_tol = 0.0;          // distance tolerance for atom movement between two time steps
    n_writefile = 1;             // number of time steps between writing the output files
    verbose_mode = "verbose";    // mute, silent, verbose
    surface_cleaner = "faces";   // method to clean surface; voronois, faces or none
    use_histclean = false;       // use histogram cleaner to get rid of sharp peaks in the solution
    cluster_anal = true;         // enable cluster analysis
    clear_output = true;         // clear output folder

    heating_mode = "none";       // method to calculate current density and temperature; none, stationary or transient
    transient_time = 0.05e-12;   // time resolution in transient temperature solver [sec]
    transient_steps = 3;         // number of iterations in transient heat equation solver

    E0 = 0;                      // long range electric field
    neumann = 0;                 // neumann boundary contition value
    message = "";                // message from the host code

    smooth_steps = 0;            // number of surface mesh smoothing iterations
    smooth_lambda = 0.5;         // lambda parameter in surface mesh smoother
    smooth_mu = -0.53;           // mu parameter in surface mesh smoother
    smooth_algorithm = "fujiwara"; // surface mesh smoother algorithm; none, laplace or fujiwara
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

    // Modify the parameters that are specified in input script
    read_command("transient_steps", transient_steps);
    read_command("transient_time", transient_time);
    read_command("t_ambient", t_ambient);
    read_command("heating_mode", heating_mode);
    read_command("smooth_steps", smooth_steps);
    read_command("smooth_lambda", smooth_lambda);
    read_command("smooth_mu", smooth_mu);
    read_command("smooth_algorithm", smooth_algorithm);
    read_command("charge_tolerance", charge_tolerance);
    read_command("phi_error", phi_error);
    read_command("n_phi", n_phi);
    read_command("surface_cleaner", surface_cleaner);
    read_command("cluster_anal", cluster_anal);
    read_command("extended_atoms", extended_atoms);
    read_command("clear_output", clear_output);
    read_command("infile", atom_file);
    read_command("latconst", latconst);
    read_command("coord_cutoff", coordination_cutoff);
    read_command("cluster_cutoff", cluster_cutoff);
    read_command("surface_thichness", surface_thichness);
    read_command("nnn", nnn);
    read_command("mesh_quality", mesh_quality);
    read_command("element_volume", element_volume);
    read_command("radius", radius);
    read_command("smooth_factor", surface_smooth_factor);
    read_command("charge_smooth_factor", charge_smooth_factor);
    read_command("refine_apex", refine_apex);
    read_command("distance_tol", distance_tol);
    read_command("box_width", box_width);
    read_command("box_height", box_height);
    read_command("bulk_height", bulk_height);    
    read_command("femocs_verbose_mode", verbose_mode);
    read_command("femocs_periodic", MODES.PERIODIC);
    read_command("write_log", MODES.WRITELOG);
    read_command("use_histclean", use_histclean);
    read_command("n_writefile", n_writefile);

    vector<double> coarse_factors = {cfactor.amplitude, (double)cfactor.r0_cylinder, (double)cfactor.r0_sphere};
    read_command("coarse_factor", coarse_factors);

    cfactor.amplitude = coarse_factors[0];
    cfactor.r0_cylinder = static_cast<int>(coarse_factors[1]);
    cfactor.r0_sphere = static_cast<int>(coarse_factors[2]);

    check_obsolete(file_name);
}

// Read the commands and their arguments from the file and store them into the buffer
void Config::parse_file(const string& file_name) {
    ifstream file(file_name);
    require(file, "File not found: " + file_name);

    string line;
    data.clear();

    // loop through the lines in a file
    while (getline(file, line)) {
        // force all the characters in a line to lower case
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        line += " "; // needed to find the end of line

        bool line_started = true;
        // store the command and its parameters from non-empty and non-pure-comment lines
        while(line.size() > 0) {
            trim(line);
            int i = line.find_first_not_of(data_symbols);
            if (i <= 0) break;

            if (line_started && line.substr(0, i) == "femocs_end") return;
            if (line_started) data.push_back({});

            data.back().push_back( line.substr(0, i) );
            line = line.substr(i);
            line_started = false;
        }
    }
}

// Check for the obsolete commands from the buffered commands
void Config::check_obsolete(const string& file_name) {
    for (vector<string> cmd : data) {
        if (obsolete_commands.find(cmd[0]) != obsolete_commands.end())
            write_verbose_msg("Command '" + cmd[0] + "' is obsolete!"
                    " You can safely remove it from " + file_name + " !");
    }
}

// Look up the parameter with string argument
int Config::read_command(string param, string& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (vector<string> str : data)
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
    for (vector<string> str : data)
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

// Look up the parameter with integer argument
int Config::read_command(string param, int& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (vector<string> str : data)
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
    for (vector<string> str : data)
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
    unsigned n_read_args = 0;
    for (vector<string> str : data)
        if (str.size() >= 2 && str[0] == param)
            for (unsigned i = 0; i < args.size() && i < (str.size()-1); ++i) {
                istringstream is(str[i+1]);
                double result;
                if (is >> result) { args[i] = result; n_read_args++; }
            }
    return n_read_args == args.size();
}

// Print the stored commands and parameters
void Config::print_data() {
    if (!MODES.VERBOSE) return;
    const int cmd_len = 20;

    for (vector<string> line : data) {
        for (string ln : line) {
            int str_len = ln.length();
            int whitespace_len = 1;
            if (cmd_len > str_len)
                whitespace_len = cmd_len - str_len;

            cout << ln << string(whitespace_len, ' ');
        }

        cout << endl;
    }
}

} // namespace femocs
