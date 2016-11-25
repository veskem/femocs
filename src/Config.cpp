/*
 * Config.cpp
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#include "Config.h"
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;
namespace femocs {

// Config constructor
Config::Config() {
    init_values();
}

// Initialize configuration parameters
const void Config::init_values() {
//    infile = "input/rough111.ckx";
//    infile = "input/mushroom2.ckx";
//    infile = "input/tower_hl2p5.ckx";
    infile = "input/nanotip_hr5.ckx";
    latconst = 2.0;

    // conf.infile = home + "input/nanotip_medium.xyz";
    // conf.latconst = 3.61;

    coord_cutoff = 3.1;         // coordination analysis cut-off radius

    nnn = 12;                    // number of nearest neighbours in bulk
    mesh_quality = "2.0";        // minimum mesh quality Tetgen is allowed to make
    nt = 4;                      // number of OpenMP threads
    radius = 14.0;               // inner radius of coarsening cylinder
    smooth_factor = 0.5;         // surface smoothing factor; bigger number gives smoother surface
    n_bins = 20;                 // number of bins in histogram smoother
    postprocess_marking = false; // make extra effort to mark correctly the vacuum nodes in shadow area
    rmin_rectancularize = latconst / 1.0; // 1.5+ for <110> simubox, 1.0 for all others

    refine_apex = false;         // refine nanotip apex
    distance_tol = 0.5*latconst;

    // Electric field is applied 100 lattice constants above the highest point of surface
    // and bulk is extended 20 lattice constants below the minimum point of surface
    zbox_above = 100 * latconst;
    zbox_below = 20 * latconst;

    cfactor.amplitude = 0.4;       // coarsening factor
    cfactor.r0_cylinder = 1.0;  // minimum distance between atoms in nanotip below apex
    cfactor.r0_sphere = 0.0;    // minimum distance between atoms in nanotip apex
}

// Remove the noise from the beginning of the string
const void Config::trim(string& str) {
    str.erase(0, str.find_first_of(comment_symbols + data_symbols));
}

// Read the configuration parameters from input script
const void Config::read_all(const string& file_name) {
    if(file_name == "") return;

    // Store the commands and their arguments
    parse_file(file_name);

    // Modify the parameters that are correctly specified in input script
    read_command("infile", infile);
    read_command("latconst", latconst);
    read_command("coord_cutoff", coord_cutoff);
    read_command("nnn", nnn);
    read_command("mesh_quality", mesh_quality);
    read_command("radius", radius);
    read_command("smooth_factor", smooth_factor);
    read_command("n_bins", n_bins);
    read_command("postprocess_marking", postprocess_marking);
    read_command("refine_apex", refine_apex);
    read_command("distance_tol", distance_tol);
    read_command("rmin_rectancularize", rmin_rectancularize);
    read_command("zbox_above", zbox_above);
    read_command("zbox_below", zbox_below);
    read_command("femocs_verbose", MODES.VERBOSE);
    read_command("femocs_writefile", MODES.WRITEFILE);

    vector<double> coarse_factors = {cfactor.amplitude, cfactor.r0_cylinder, cfactor.r0_sphere};
    read_command("coarse_factor", coarse_factors);

    cfactor.amplitude = coarse_factors[0];
    cfactor.r0_cylinder = coarse_factors[1];
    cfactor.r0_sphere = coarse_factors[2];
}

const void Config::parse_file(const string& file_name) {
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

            if (line_started && line.substr(0, i) == "end") return;
            if (line_started) data.push_back({});

            data.back().push_back( line.substr(0, i) );
            line = line.substr(i);
            line_started = false;
        }
    }
}

// Look up the parameter with string argument
const int Config::read_command(const string& param, string& arg) {
    // loop through all the commands that were found from input script
    for (vector<string> str : data)
        if (str.size() >= 2 && str[0] == param) {
            arg = str[1]; return 0;
        }
    return 1;
}

// Look up the parameter with boolean argument
const int Config::read_command(const string& param, bool& arg) {
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
const int Config::read_command(const string& param, int& arg) {
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
const int Config::read_command(const string& param, double& arg) {
    // loop through all the commands that were found from input script
    for (vector<string> str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); double result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

const int Config::read_command(const string& param, vector<double>& args) {
    // loop through all the commands that were found from input script
    int n_read_args = 0;
    for (vector<string> str : data)
        if (str.size() >= 2 && str[0] == param)
            for (int i = 0; i < args.size() && i < (str.size()-1); ++i) {
                istringstream is(str[i+1]);
                double result;
                if (is >> result) { args[i] = result; n_read_args++; }
            }
    return n_read_args == args.size();
}
// Print the stored commands and parameters
const void Config::print_data() {
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
