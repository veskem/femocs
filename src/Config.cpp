/*
 * Config.cpp
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#include "Config.h"
#include <fstream>
#include <sstream>

using namespace std;
namespace femocs {

Config::Config() {
    init_values();
}

const void Config::init_values() {
//    infile = "input/rough111.ckx";
//    infile = "input/mushroom2.ckx";
//    infile = "input/tower_hl2p5.ckx";
    infile = "input/nanotip_hr5.ckx";
    latconst = 2.0;

    // conf.infile = home + "input/nanotip_medium.xyz";
    // conf.latconst = 3.61;

    coord_cutoff = 3.1;         // coordination analysis cut-off radius

    nnn = 12;                   // number of nearest neighbours in bulk
    mesh_quality = "2.0";       // minimum mesh quality Tetgen is allowed to make
    nt = 4;                     // number of OpenMP threads
    radius = 10.0;              // inner radius of coarsening cylinder
    coarse_factor = 0.4;        // coarsening factor; bigger number gives coarser surface
    smooth_factor = 0.5;        // surface smoothing factor; bigger number gives smoother surface
    n_bins = 20;                // number of bins in histogram smoother
    postprocess_marking = true; // make extra effort to mark correctly the vacuum nodes in shadow area
    rmin_rectancularize = latconst / 1.0; // 1.5+ for <110> simubox, 1.0 for all others

    refine_apex = false;        // refine nanotip apex
    significant_distance = 0.5*latconst;

    // Electric field is applied 100 lattice constants above the highest point of surface
    // and bulk is extended 20 lattice constants below the minimum point of surface
    zbox_above = 100 * latconst;
    zbox_below = 20 * latconst;
}

const void Config::trim(string& str) {
    str.erase(0, str.find_first_not_of(" =\t"));
    str.erase(str.find_last_not_of(" =\t")+1);
}

const void Config::read_all(const string& file_name) {
    ifstream file(file_name);
    string line;
    vector<vector<string>> data;

    while (getline(file, line)) {
        data.push_back({});
        while(line.size() > 0) {
            trim(line);
            int i = line.find_first_of(" #=\t");
            if (i <= 0) break;
            data[data.size()-1].push_back( line.substr(0, i) );
            line = line.substr(i);
        }
    }

//    for (int i = 0; i < data.size(); ++i) {
//        for (int j = 0; j < data[i].size(); ++j)
//            cout << data[i][j] << "\t";
//        cout << endl;
//    }
//    istringstream ss( line.substr(line.find("=") + 1) );
//    cout << ss;
//    if (line.find("num") != -1)
//        ss >> config.num;
//    else if (line.find("str") != -1)
//        ss >> config.str;
//    else if (line.find("flt") != -1)
//        ss >> config.flt;
}

const void Config::read_parameter(const string& param, string& value) {

}
const void Config::read_parameter(const string& param, bool& value) {

}

const void Config::read_parameter(const string& param, int& value) {

}

const void Config::read_parameter(const string& param, double& value) {

}

} // namespace femocs

