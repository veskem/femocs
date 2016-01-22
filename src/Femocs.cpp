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
#include "Mesher.h"
#include "SurfaceExtractor.h"
#include "Vacuum.h"

using namespace std;
namespace femocs {

// Femocs constructor
Femocs::Femocs(const string& file_name) {
    conf = parse_input_script(file_name);
    simucell = init_simucell();
}

// Get simulation parameters
const Femocs::Config Femocs::parse_input_script(const string& fileName) const {
    Config conf;
    conf.latconst = 1.0; //3.51;       // lattice constant
    conf.coord_cutoff = 3.3;	// coordination analysis cut off radius
    conf.tetgen_cutoff = 14.2; //8.0;			// max_length^2 of tetrahedra's edge
    conf.nnn = 12;						    // number of nearest neighbours in bulk
    conf.extracter = "coordination";	    // surface extraction algorithm
    conf.mesher = "tetgen";				    // mesher algorithm
    conf.nt = 4;						    // number of OpenMP threads

    return conf;
}

// Initialise general data about simulation cell
const Femocs::SimuCell Femocs::init_simucell() const {
    SimuCell sc;
    sc.xmin = sc.xmax = sc.ymin = sc.ymax = sc.zmin = sc.zmax = 0;
    sc.zmaxbox = 100.0;
    return sc;
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

// ===============================================
// ========== Start of Femocs interface ==========
// ***********************************************
int main(int argc, char* argv[]) {
    using namespace femocs;
//    string in_file_name = "input/nanotip_medium.xyz";
    string in_file_name = "input/mushroomish_surface.ckx";
    string surface_file = "output/surface.xyz";
    string bulk_file = "output/bulk.xyz";
    string vacuum_file = "output/vacuum.xyz";

    int nt;
    double t0;
    Femocs femocs("");

    // If input file specified on command line, use that instead of default
    if (argc == 2) in_file_name = argv[1];
    // The same with number of OpenMP threads
    nt = femocs.conf.nt;
    if (argc == 3) {
        in_file_name = argv[1];
        nt = stod(argv[2]);
    }

    t0 = start_msg("=== Reading atoms from " + in_file_name);
    AtomReader reader;
    reader.import_file(in_file_name, &femocs.simucell);
    end_msg(t0);

    t0 = start_msg("=== Extracting surface...");
    SurfaceExtractor surf_extractor(femocs.conf.extracter, femocs.conf.coord_cutoff,
            femocs.conf.latconst, femocs.conf.nnn, nt);
    auto surf = surf_extractor.extract_surface(&reader.data, &femocs.simucell);
    end_msg(t0);

    t0 = start_msg("=== Extracting bulk...");
    SurfaceExtractor bulk_extractor(femocs.conf.extracter, 0, femocs.conf.latconst, femocs.conf.nnn,
            nt);
    auto bulk = bulk_extractor.extract_surface(&reader.data, &femocs.simucell);
    end_msg(t0);

    t0 = start_msg("=== Generating vacuum...");
    Vacuum vacuum;
    vacuum.generate(&femocs.simucell, surf);
    end_msg(t0);

//    t0 = start_msg("=== Writing surface, bulk and vacuum to output/...");
//    surf->output(surface_file);
//    bulk->output(bulk_file);
//    vacuum.output(vacuum_file);
//    end_msg(t0);

    t0 = start_msg("=== Generating surface mesh...\n");
    Mesher mesher(femocs.conf.mesher);
    mesher.generate_mesh(&femocs.simucell, bulk, surf, &vacuum, "Q");
    end_msg(t0);

//    t0 = start_msg("=== Cleaning surface mesh...\n");
//    mesher.clean_mesh("r", femocs.conf.tetgen_cutoff);
//    end_msg(t0);

//    t0 = start_msg("=== Refining mesh...");
//    mesher.refine("rq1.414a100Y");
//    end_msg(t0);

    t0 = start_msg("=== Marking mesh...");
    mesher.mark_mesh(&femocs.simucell, "r");
    end_msg(t0);

    t0 = start_msg("=== Outputting mesh...");
    mesher.write_faces("output/faces.vtk");
    mesher.write_elems("output/elems.vtk");
    end_msg(t0);

    cout << "\n======= Femocs finished! =======\n";

    return 0;
}
