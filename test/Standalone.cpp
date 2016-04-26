/*
 * Standalone.cpp
 *
 *  Created on: 19.04.2016
 *      Author: veske
 */

#include "Femocs.h"

using namespace std;

int main(int argc, char* argv[]) {
    Femocs femocs("/path/to/input.script");

    // If input file specified on command line, use that instead of default
    if (argc >= 2) femocs.conf.infile = argv[1];
    // The same with number of OpenMP threads
    if (argc >= 3) femocs.conf.nt = stod(argv[2]);

    // Make some dummy variables that resemble Helmod input/output format
    const int N = 2;
    int i, j;
    double E0 = 0.0;
    double x[2] = {0.0, 0.0};

    double grid_spacing[3] = {femocs.conf.latconst, femocs.conf.latconst, femocs.conf.latconst};
    int*** BC = new int**[N];
    double*** phi = new double**[N];

    for (i = 0; i < N; ++i) {
        BC[i] = new int*[N];
        phi[i] = new double*[N];
        for (j = 0; j < N; ++j) {
            BC[i][j] = new int[N];
            phi[i][j] = new double[N];
        }
    }

    // Run the actual script
//    femocs.run(E0, BC, phi_guess, grid_spacing);

    femocs.import_atoms(0, x, x, x);
    femocs.run(E0, phi);

    return 0;
}
