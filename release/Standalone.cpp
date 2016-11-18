/*
 * Standalone.cpp
 *
 *  Created on: 19.04.2016
 *      Author: veske
 */

#include "Femocs.h"
#include "Macros.h"

using namespace std;

//int main(int argc, char* argv[]) {
int main() {
    cout << "Standalone started!" << endl;

    int n_atoms = 46053;
    int n_points = 100;
    double x[n_points];
    double y[n_points];
    double z[n_points];
    double phi[n_points];
    int flag[n_points];

    double Ex[n_atoms];
    double Ey[n_atoms];
    double Ez[n_atoms];
    double Enorm[n_atoms];

    for (int i = 0; i < n_points; ++i) {
        x[i] = 0; y[i] = 0; z[i] = 1.0 * i;
    }

    int success = 0;

    femocs::Femocs femocs("input/md.in");
    success += femocs.import_atoms("");
    success += femocs.run(0.1, "");
    success += femocs.export_elfield(n_atoms, Ex, Ey, Ez, Enorm);
    success += femocs.interpolate_phi(n_points, x, y, z, phi, flag);
    success += femocs.interpolate_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);

    cout << "Standalone result: " << success << endl;

    return 0;
}
