/*
 * Standalone.cpp
 *
 *  Created on: 19.04.2016
 *      Author: veske
 */

#include "Femocs.h"

using namespace std;

//int main(int argc, char* argv[]) {
int main() {
    Femocs femocs("/path/to/input.script");

//    // If input file specified on command line, use that instead of default
//    if (argc >= 2) femocs.conf.infile = argv[1];
//    // The same with number of OpenMP threads
//    if (argc >= 3) femocs.conf.nt = stod(argv[2]);

    // Make some dummy variables that resemble Helmod input/output format
    double E0 = 0.1;
    double x[2] = {0.0, 0.0};
    int type[2] = {1, 1};

    // Run the actual script
    femocs.import_atoms(0, x, x, x, type);
    femocs.run2(E0);

    return 0;
}
