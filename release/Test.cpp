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

using namespace std;

void print_progress(const string& message, const bool contition) {
    cout << message << ":\t";
    if (contition) cout << "passed" << endl;
    else cout << "failed" << endl;
}

void write_defaults(ofstream &file) {
    file << "heating = false"            << endl;
    file << "clear_output = true"        << endl;
    file << "smooth_factor = 0.0"        << endl;
    file << "charge_smooth_factor = 1.0" << endl;
    file << "postprocess_marking = true" << endl;
    file << "distance_tol = 0.16"        << endl;
    file << "n_writefile = 1"            << endl;
    file << "femocs_verbose = true"      << endl;
    file << "coord_cutoff = 3.1"         << endl;
    file << "latconst = 3.61"            << endl;
    file << "use_histclean = false"      << endl;
    file << "surface_cleaner = voronois" << endl;
}

void write_hr5(ofstream &file) {
    file << "infile = input/nanotip_hr5.xyz" << endl;
    file << "coarse_factor = 0.3 4 2"        << endl;
    file << "radius = 14.0"                  << endl;
}

void write_rectangle(ofstream &file) {
    file << "infile = input/nanotip_rectangle.xyz" << endl;
    file << "coarse_factor = 0.3 4 2"        << endl;
    file << "radius = 14.0"                  << endl;
}

void write_mdsmall(ofstream &file) {
    file << "infile = input/nanotip_small.xyz" << endl;
    file << "coarse_factor = 0.3 4 2"    << endl;
    file << "radius = 16.0"              << endl;
    file << "box_width = 4.0"            << endl;
    file << "box_height = 3.5"           << endl;
}

void write_mdbig(ofstream &file) {
    file << "infile = input/nanotip_big.xyz" << endl;
    file << "coarse_factor = 0.3 4 2" << endl;
    file << "radius = 16.0"           << endl;
}

void write_kmcsmall(ofstream &file) {
    file << "infile = input/mushroom1.ckx" << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "latconst = 2.0"          << endl;
    file << "radius = 11.0"           << endl;
}

void write_kmcbig(ofstream &file) {
    file << "infile = input/mushroom2.ckx" << endl;
    file << "coarse_factor = 0.4 6 4" << endl;
    file << "latconst = 2.0"          << endl;
    file << "radius = 20.0"           << endl;
}

void write_stretch(ofstream &file) {
    file << "infile = input/nanotip_big.xyz" << endl;
    file << "coarse_factor = 0.3 4 2" << endl;
    file << "radius = 16.0"           << endl;
    file << "box_width = 4.0"         << endl;
    file << "box_height = 3.5"        << endl;
    file << "bulk_height = 20"        << endl;
}

void write_extend(ofstream &file) {
    file << "extended_atoms = input/extension.xyz" << endl;
    file << "infile = input/apex.xyz"              << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "femocs_periodic = false" << endl;
    file << "radius = 70.0"           << endl;
}

void write_cluster(ofstream &file) {
    file << "infile = input/clusters.xyz" << endl;
    file << "coarse_factor = 0.3 6 4" << endl;
    file << "radius = 12.0"           << endl;
}

void write_molten(ofstream &file) {
    file << "infile = input/nanotip_molten.xyz" << endl;
    file << "coarse_factor = 0.3 6 4"  << endl;
    file << "radius = 45.0"            << endl;
    file << "surface_thichness = 4.65" << endl;
}

void write_moltenbig(ofstream &file) {
    file << "infile = input/nanotip_molten.ckx" << endl;
    file << "coarse_factor = 0.3 6 4"  << endl;
    file << "radius = 65.0"            << endl;
    file << "surface_thichness = 4.65" << endl;
}

int main(int argc, char **argv) {
    string filename = "input/md.in";
    string mode = "default";
    char arg[128];

    if (argc >= 2) {
        filename = "tmpfile";

        sscanf(argv[1], "%s", arg);
        mode = string(arg);

        ofstream file(filename);

        if      (mode == "kmcsmall")  write_kmcsmall(file);
        else if (mode == "kmcbig")    write_kmcbig(file);
        else if (mode == "hr5")       write_hr5(file);
        else if (mode == "rectangle") write_rectangle(file);
        else if (mode == "mdsmall")   write_mdsmall(file);
        else if (mode == "mdbig")     write_mdbig(file);
        else if (mode == "stretch")   write_stretch(file);
        else if (mode == "extend")    write_extend(file);
        else if (mode == "cluster")   write_cluster(file);
        else if (mode == "molten")    write_molten(file);
        else if (mode == "moltenbig") write_moltenbig(file);

        else {
            printf("Usage:\n");
            printf("  no-arg      configuration is obtained from input/md.in\n");
            printf("  kmcsmall    small kMC nanotip\n");
            printf("  kmcbig      big kMC nanotip\n");
            printf("  mdsmall     small MD nanotip\n");
            printf("  mdbig       big MD nanotip\n");
            printf("  hr5         symmetric nanotip with aspect ratio 5\n");
            printf("  rectangle   symmetric nanotip with rectangular substrate\n");
            printf("  stretch     stretch the substrate of small MD nanotip\n");
            printf("  extend      extend the system below the MD apex\n");
            printf("  cluster     MD nanotip with clusters\n");
            printf("  molten      nanotip with molten apex on top of thin rod\n");
            printf("  moltenbig   symmetric MD nanotip with molten apex\n");

            file.close();
            system("rm -rf tmpfile");
            exit(0);
        }

        write_defaults(file);
        file.close();
    }

    int n_atoms = 1000;
    int n_points = 100;
    int flag[n_points];
    double x[n_points], y[n_points], z[n_points], phi[n_points];
    double Ex[n_atoms], Ey[n_atoms], Ez[n_atoms], Enorm[n_atoms], T[n_atoms], xq[4*n_atoms];

    for (int i = 0; i < n_points; ++i) {
        x[i] = 0; y[i] = 0; z[i] = 1.0 * i;
    }

    cout << "\nRunning FEMOCS test program in a mode   " << mode << endl;

    int success = 0;
    femocs::Femocs femocs(filename);

    system("rm -rf tmpfile");

    success += femocs.import_atoms("");
    success += femocs.run(-0.1, "");
    success += femocs.export_elfield(n_atoms, Ex, Ey, Ez, Enorm);
    success += femocs.export_temperature(n_atoms, T);
    success += femocs.export_charge_and_force(n_atoms, xq);
    success += femocs.interpolate_phi(n_points, x, y, z, phi, flag);
    success += femocs.interpolate_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);

    cout << endl;
    print_progress("full run of Femocs", success==0);

    string cmd1 = "Smooth_Factor"; double arg1 = -1.0;
    success = femocs.parse_command(cmd1, &arg1);
    print_progress("reading "+cmd1, arg1 != -1.0);

    return 0;
}
