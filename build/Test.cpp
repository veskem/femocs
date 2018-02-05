/*
 * Test program to demonstrate and test FEMOCS with various geometries.
 * Testing mode is picked from command line argument -
 * without argument the configuration script from input folder is used,
 * but with cmd line parameter other pre-tuned modes can be chosen.
 *
 *  Created on: 19.04.2016
 *      Author: veske
 */

#include "Femocs.h"
#include "Macros.h"

#include <stdlib.h>

using namespace std;

void print_progress(const string& message, const bool contition) {
    cout << message << ":  ";
    if (contition) cout << "passed" << endl;
    else cout << "failed" << endl;
}

void write_defaults(ofstream &file) {
    file << "mesh_quality = 1.8"         << endl;
    file << "heating_mode = none"        << endl;
    file << "write_log = true"           << endl;
    file << "clear_output = true"        << endl;
    file << "surface_smooth_factor= 0.1" << endl;
    file << "charge_smooth_factor = 1.0" << endl;
    file << "distance_tol = 0.0"         << endl;
    file << "n_writefile = 1"            << endl;
    file << "use_histclean = false"      << endl;
    file << "use_rdf = false"            << endl;
    file << "clean_surface = true"    << endl;
    file << "surface_thickness = 3.1"    << endl;
    file << "coord_cutoff = 3.1"         << endl;
    file << "latconst = 3.61"            << endl;
    file << "femocs_verbose_mode = verbose" << endl;
    file << "smooth_steps = 3"           << endl;
    file << "smooth_algorithm = laplace" << endl;
    file << "elfield = -0.07"            << endl;
    file << "interpolation_rank = 1"     << endl;
}

void write_rectangle(ofstream &file) {
    file << "infile = in/nanotip_rectangle.xyz" << endl;
    file << "coarse_factor = 0.3 4 2"        << endl;
    file << "radius = 14.0"                  << endl;
}

void write_mdsmall(ofstream &file) {
    file << "infile = in/nanotip_small.xyz" << endl;
    file << "coarse_factor = 0.3 4 2"    << endl;
    file << "radius = 16.0"              << endl;
    file << "box_width = 4.0"            << endl;
    file << "box_height = 3.5"           << endl;
}

void write_wobble(ofstream &file) {
    write_mdsmall(file);
}

void write_mdbig(ofstream &file) {
    file << "infile = in/nanotip_big.xyz" << endl;
    file << "coarse_factor = 0.3 4 2" << endl;
    file << "radius = 16.0"           << endl;
}

void write_kmcsmall(ofstream &file) {
    file << "infile = in/mushroom1.ckx" << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "latconst = 2.0"          << endl;
    file << "radius = 11.0"           << endl;
    file << "clean_surface = false"  << endl;
}

void write_kmcbig(ofstream &file) {
    file << "infile = in/mushroom2.ckx" << endl;
    file << "coarse_factor = 0.4 6 4" << endl;
    file << "latconst = 2.0"          << endl;
    file << "radius = 20.0"           << endl;
    file << "clean_surface = false"  << endl;
    file << "mesh_quality = 1.6"      << endl;
}

void write_kmcregular(ofstream &file) {
    file << "infile = in/kmc_regular.ckx" << endl;
    file << "coarse_factor = 0.4 6 4" << endl;
    file << "latconst = 3.6935"       << endl;
    file << "radius = 42.0"           << endl;
    file << "box_width = 5.0"         << endl;
    file << "box_height = 5.0"        << endl;
    file << "clean_surface = false"  << endl;
}

void write_stretch(ofstream &file) {
    file << "infile = in/nanotip_big.xyz" << endl;
    file << "coarse_factor = 0.3 4 2" << endl;
    file << "radius = 16.0"           << endl;
    file << "box_width = 4.0"         << endl;
    file << "box_height = 3.5"        << endl;
    file << "bulk_height = 20"        << endl;
}

void write_extend(ofstream &file) {
    file << "extended_atoms = in/extension.xyz" << endl;
    file << "infile = in/apex.ckx"              << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "femocs_periodic = false" << endl;
    file << "radius = 70.0"           << endl;
}

void write_tablet(ofstream &file) {
    file << "extended_atoms = in/extension.xyz" << endl;
    file << "infile = in/tablet.ckx"              << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "femocs_periodic = false" << endl;
    file << "radius = 70.0"           << endl;
}

void write_cluster(ofstream &file) {
    file << "infile = in/clusters.xyz" << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "radius = 12.0"           << endl;
}

void write_molten(ofstream &file) {
    file << "infile = in/nanotip_molten.xyz" << endl;
    file << "coarse_factor = 0.3 6 4"  << endl;
    file << "radius = 65.0"            << endl;
    file << "box_width = 5.0"         << endl;
}

void write_moltenbig(ofstream &file) {
    file << "infile = in/nanotip_molten.ckx" << endl;
    file << "coarse_factor = 0.4 8 3"  << endl;
    file << "radius = 45.0"            << endl;
}

void write_tip100(ofstream &file) {
    file << "infile = in/tip100.ckx"   << endl;
    file << "coarse_factor = 0.4 8 3"  << endl;
    file << "radius = 45.0"            << endl;
}

void write_tip110(ofstream &file) {
    file << "infile = in/tip110.ckx"   << endl;
    file << "coarse_factor = 0.4 8 3"  << endl;
    file << "radius = 45.0"            << endl;
}

void write_tip111(ofstream &file) {
    file << "infile = in/tip111.ckx"   << endl;
    file << "coarse_factor = 0.4 8 3"  << endl;
    file << "radius = 45.0"            << endl;
}

void write_generate(ofstream &file) {
    file << "infile = in/nanotip_small.xyz" << endl;
    file << "coarse_factor = 0.2 1 1"  << endl;
    file << "radius = 31.0"            << endl;
    file << "box_width = 10.0"         << endl;
    file << "box_height = 10.0"        << endl;
    file << "clean_surface = false"   << endl;
}

void read_xyz(const string &file_name, double* x, double* y, double* z) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    int n_atoms, type, id;
    string elem, line;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;

    getline(in_file, line);     // Skip comments line

    id = -1;
    // keep storing values from the text file as long as data exists:
    while (++id < n_atoms && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> elem >> x[id] >> y[id] >> z[id] >> type;
    }
}

void read_ckx(const string &file_name, double* x, double* y, double* z) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    int n_atoms, type, id;
    string line;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;

    getline(in_file, line);    // Skip comments line

    id = -1;
    // keep storing values from the text file as long as data exists:
    while (++id < n_atoms && getline(in_file, line)) {
        iss.clear();
        iss.str(line);
        iss >> type >> x[id] >> y[id] >> z[id];
    }
}

void read_atoms(const string& file_name, double* x, double* y, double* z) {
    string file_type = get_file_type(file_name);

    if (file_type == "xyz")
        read_xyz(file_name, x, y, z);
    else if (file_type == "ckx")
        read_ckx(file_name, x, y, z);
    else
        require(false, "Unsupported file type: " + file_type);
}

void read_n_atoms(const string& file_name, int& n_atoms) {
    ifstream in_file(file_name, ios::in);
    require(in_file.is_open(), "Did not find a file " + file_name);

    string line;
    istringstream iss;

    getline(in_file, line);     // Read number of atoms
    iss.clear();
    iss.str(line);
    iss >> n_atoms;
}

int main(int argc, char **argv) {
    string filename = "in/md.in";
    string mode = "default";
    char arg[128];
    int success = 0;

    // read the running mode
    if (argc >= 2) {
        filename = "tmpfile";

        sscanf(argv[1], "%s", arg);
        mode = string(arg);

        ofstream file(filename);

        if      (mode == "kmcsmall")  write_kmcsmall(file);
        else if (mode == "kmcbig")    write_kmcbig(file);
        else if (mode == "kmcregular") write_kmcregular(file);
        else if (mode == "tip100")    write_tip100(file);
        else if (mode == "tip110")    write_tip110(file);
        else if (mode == "tip111")    write_tip111(file);
        else if (mode == "rectangle") write_rectangle(file);
        else if (mode == "mdsmall")   write_mdsmall(file);
        else if (mode == "mdbig")     write_mdbig(file);
        else if (mode == "stretch")   write_stretch(file);
        else if (mode == "extend")    write_extend(file);
        else if (mode == "tablet")    write_tablet(file);
        else if (mode == "cluster")   write_cluster(file);
        else if (mode == "molten")    write_molten(file);
        else if (mode == "moltenbig") write_moltenbig(file);
        else if (mode == "generate")  write_generate(file);
        else if (mode == "wobble")    write_wobble(file);

        else {
            printf("Usage:\n");
            printf("  no-arg      configuration is obtained from in/md.in\n");
            printf("  kmcsmall    small kMC nanotip\n");
            printf("  kmcbig      big kMC nanotip\n");
            printf("  kmcregular  symmetric kMC nanotip\n");
            printf("  mdsmall     small MD nanotip\n");
            printf("  mdbig       big MD nanotip\n");
            printf("  tip100      symmetric nanotip with h/r = 5 and <100> orientation\n");
            printf("  tip110      symmetric nanotip with h/r = 5 and <110> orientation\n");
            printf("  tip111      symmetric nanotip with h/r = 5 and <111> orientation\n");
            printf("  rectangle   symmetric nanotip with rectangular substrate\n");
            printf("  stretch     stretch the substrate of small MD nanotip\n");
            printf("  extend      extend the system below the round MD apex\n");
            printf("  tablet      extend the system below the tablet shaped MD apex\n");
            printf("  cluster     MD nanotip with clusters\n");
            printf("  molten      nanotip with molten apex on top of thin rod\n");
            printf("  moltenbig   symmetric MD nanotip with molten apex\n");
            printf("  generate    generate and use perfectly symmetric nanotip without crystallographic properties\n");
            printf("  wobble      read small MD nanotip and add random noise to emulate real MD simulation\n");

            file.close();
            success = system("rm -rf tmpfile");
            exit(0);
        }

        write_defaults(file);
        file.close();
    }

    cout << "\n> running FEMOCS test program in a mode:  " << mode << endl;

    femocs::Femocs femocs(filename);
    success = system("rm -rf tmpfile");

    string cmd1 = "infile"; string infile = "";
    success = femocs.parse_command(cmd1, infile);
    print_progress("\n> reading " + cmd1, infile != "");

    // determine number of iterations
    int n_iterations = 1;
    if (argc >= 3) {
        n_iterations = atoi(argv[2]);
    }

    int n_atoms = 0, n_points = 100;
    read_n_atoms(infile, n_atoms);

    int* flag  = (int*)    malloc(n_atoms * sizeof(int));
    double* x  = (double*) malloc(n_atoms * sizeof(double));
    double* y  = (double*) malloc(n_atoms * sizeof(double));
    double* z  = (double*) malloc(n_atoms * sizeof(double));
    double* phi= (double*) malloc(n_atoms * sizeof(double));
    double* Ex = (double*) malloc(n_atoms * sizeof(double));
    double* Ey = (double*) malloc(n_atoms * sizeof(double));
    double* Ez = (double*) malloc(n_atoms * sizeof(double));
    double* En = (double*) malloc(n_atoms * sizeof(double));
    double* T  = (double*) malloc(n_atoms * sizeof(double));
    double* xq = (double*) malloc(4 * n_atoms * sizeof(double));

    read_atoms(infile, x, y, z);

    for (int i = 1; i <= n_iterations; ++i) {
        if (n_iterations > 1) cout << "\n> iteration " << i << endl;

        if (mode == "generate") {
            if (argc >= 3)
                success += femocs.generate_nanotip(0, 30, atof(argv[2]));
            else
                success += femocs.generate_nanotip(0, 30);
        } else if (mode == "wobble")
            success += femocs.import_atoms(infile, true);
        else
            success += femocs.import_atoms(infile);

        success += femocs.run();
        success += femocs.export_elfield(0, Ex, Ey, Ez, En);
//        success += femocs.export_temperature(n_atoms, T);
        success += femocs.export_charge_and_force(n_atoms, xq);
//        success += femocs.interpolate_elfield(n_atoms, x, y, z, Ex, Ey, Ez, En, flag);
//        success += femocs.interpolate_phi(n_atoms, x, y, z, phi, flag);
    }

    print_progress("\n> full run of Femocs", success == 0);

    free(flag);
    free(x); free(y); free(z);
    free(Ex); free(Ey); free(Ez); free(En);
    free(phi); free(T); free(xq);

    return 0;
}
