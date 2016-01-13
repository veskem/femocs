/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"

#include <omp.h>
#include <iostream>
#include <memory>

#include "AtomReader.h"
#include "Medium.h"
#include "Mesher.h"
#include "SurfaceExtractor.h"

using namespace std;
namespace femocs {

// Femocs constructor
Femocs::Femocs(const string& file_name) :
        conf(parse_input_script(file_name)) {};

// Get simulation parameters
const Femocs::Config Femocs::parse_input_script(const string& fileName) const {
    Config conf;
    conf.coord_cutoff = 3.5; // 0.0;			// coordination analysis cut off radius
    conf.tetgen_cutoff = 80.0; //8.0;			// max_length^2 of tetrahedra's edge
    conf.nnn = 12;						// number of nearest neighbours in bulk
    conf.extracter = "coordination";	// surface extraction algorithm
    conf.mesher = "tetgen";				// mesher algorithm
    conf.smoother = "laplace";			// mesh solver algorithm
    conf.nt = 1;						// number of OpenMP threads

    return conf;
}

} /* namespace femocs */

// Inliners to print informative and timing messages about processes
inline double start_msg(string s) {
    cout << endl << s;
    cout.flush();
    return omp_get_wtime();
}
inline void end_msg(double t0) {
    cout << ", time: " << omp_get_wtime() - t0 << endl;
}

// Main function of Femocs code
int main(int argc, char* argv[]) {
    using namespace femocs;
    string inFileName = "input/nanotip_medium.xyz";
    string outFileName = "output/surface.xyz";

    int nt;
    double t0;
    Femocs femocs("");

    // If input file specified on command line, use that instead of default
    if (argc == 2) inFileName = argv[1];
    // The same with number of OpenMP threads
    nt = femocs.conf.nt;
    if (argc == 3) {
        inFileName = argv[1];
        nt = stod(argv[2]);
    }

    t0 = start_msg("=== Reading atoms from " + inFileName);
    AtomReader reader(inFileName);
    end_msg(t0);

    t0 = start_msg("=== Extracting surface...");
    SurfaceExtractor extractor(femocs.conf.extracter, femocs.conf.coord_cutoff, femocs.conf.nnn,
            nt);
    auto surf = extractor.extract_surface(reader.data);
    end_msg(t0);

    t0 = start_msg("=== Writing surface to " + outFileName);
    surf->output(outFileName);
    end_msg(t0);

    t0 = start_msg("=== Generating surface mesh...\n");
    Mesher mesher(femocs.conf.mesher);
    mesher.generate_mesh(surf, "Q");
    end_msg(t0);

    t0 = start_msg("=== Cleaning surface mesh...\n");
    mesher.clean_mesh("rQ", femocs.conf.tetgen_cutoff);
    end_msg(t0);

//	t0 = startMsg("=== Smoothing surface mesh...");
//	mesher.smooth_mesh(5);
//	endMsg(t0);

    t0 = start_msg("=== Outputting surface mesh...");
    mesher.output_mesh();
    end_msg(t0);

    cout << "\n======= Femocs finished! =======\n";

    return 0;
}
