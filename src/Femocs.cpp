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

    conf.infile = "input/mushroom1.ckx";
    conf.latconst = 1.0;        // lattice constant
    conf.coord_cutoff = 3.3;    // coordination analysis cut off radius
    conf.tetgen_cutoff = 4.1;   // max_length^2 of tetrahedra's edge

//    conf.infile = "input/nanotip_medium.xyz";
//    conf.latconst = 3.51;       // lattice constant
//    conf.coord_cutoff = 3.4;	// coordination analysis cut off radius
//    conf.tetgen_cutoff = 16.0;  // max_length^2 of tetrahedra's edge

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
    sc.zmaxbox = 30*conf.latconst; // Electric field will be applied that much higer from the highest coordinate of simubox
    sc.zminbox = 0.0;
    return sc;
}

// Inliners to print messages and execution times about processes
inline double start_msg(string s) {
#if DEBUGMODE
    cout << endl << s;
    cout.flush();
#endif
    return omp_get_wtime();
}
inline void end_msg(double t0) {
#if DEBUGMODE
    cout << ", time: " << omp_get_wtime() - t0 << endl;
#endif
}

// Workhorse function to run Femocs simulation
const void Femocs::run_femocs(const double E0, double*** BC, double*** phi_guess, const double* grid_spacing) {
    AtomReader reader;
    double t0;
    // Max tetrahedron volume is ~1000x the volume of regular tetrahedron with edge == latconst
    string Vmax = to_string(  (int)(1000.0*0.118*pow(conf.latconst,3.0) ) );

#if not LIBRARYMODE
    t0 = start_msg("=== Reading atoms from " + conf.infile);
    reader.import_file(conf.infile, &simucell);
    end_msg(t0);
#elif HELMODMODE
    reader.import_helmod(E0, BC, phi_guess,...);
#elif KIMOCSMODE
    reader.import_kimocs();
#else
    cout << "Check definitions in Femocs.h !" << endl;
    exit(EXIT_FAILURE);
#endif

    t0 = start_msg("=== Extracting surface...");
    SurfaceExtractor surf_extractor(&conf);
    auto surf = surf_extractor.extract_surface(&reader.data, &simucell);
    end_msg(t0);

    t0 = start_msg("=== Extracting bulk...");
    SurfaceExtractor bulk_extractor(&conf);
//    auto bulk = bulk_extractor.extract_bulk_reduced(surf, &femocs.simucell);
    auto bulk = bulk_extractor.extract_bulk(&reader.data, &simucell);
    end_msg(t0);

    t0 = start_msg("=== Generating vacuum...");
    Vacuum vacuum;
    vacuum.generate_simple(&simucell);
    end_msg(t0);

    t0 = start_msg("=== Writing surface, bulk and vacuum to output/...");
    surf->output("output/surface.xyz");
    bulk->output("output/bulk.xyz");
    vacuum.output("output/vacuum.xyz");
    end_msg(t0);

// ===============================================

/*
 * Vahur's recipe for making tetgen mesh:
 *
 * generate mesh1 for bulk material
 * remove elements
 * rebuild mesh1
 *
 * generate mesh2 with rebuilt elements from mesh1 with -r
 * remove faces
 * rebuild mesh2
 *
 * now mesh2 = bulk elements + bulk faces
 * from mesh2, remove faces on simucell edges
 * now mesh2 = material-vacuum interface faces
 *
 * make mesh3 with material atoms + 4 points on top of simucell
 * now mesh3 = face list for the outer edges of simucell
 *
 * make mesh4 = mesh3 faces + mesh2 faces
 * mark mesh4  x,y,z min-max edges and material surface
 *
 * rebuild mesh4 with max volume and min quality restrictions
 *
 */

    t0 = start_msg("=== Step1...");

    Mesher mesher1(conf.mesher);
    auto mesh1 = mesher1.get_bulk_mesh(bulk, "Q");
    mesh1->write_faces("output/faces0.vtk");
    mesh1->write_elems("output/elems0.vtk");

    mesher1.clean_elems(mesh1, conf.tetgen_cutoff, "rQ");
    mesh1->write_faces("output/faces1.vtk");
    mesh1->write_elems("output/elems1.vtk");

    mesher1.mark_faces(mesh1, surf, &simucell);
    mesher1.calc_statistics(mesh1);
    mesh1->write_faces("output/faces2.vtk");
    mesh1->write_elems("output/elems2.vtk");

    end_msg(t0);
    t0 = start_msg("=== Step2...");

    Mesher mesher2(conf.mesher);
    auto mesh2 = mesher2.get_volume_mesh(bulk, &vacuum, "Q");
    mesh2->write_faces("output/faces3.vtk");
    mesh2->write_elems("output/elems3.vtk");

    mesh2->recalc("rq1.414a"+Vmax);
    mesh2->write_faces("output/faces4.vtk");
    mesh2->write_elems("output/elems4.vtk");

    end_msg(t0);
    t0 = start_msg("=== Step3...");

    Mesher mesher3(conf.mesher);
    auto mesh3 = mesher3.get_union_mesh(mesh1, mesh2, &simucell);
    mesh3->write_faces("output/faces5.vtk");
    mesh3->write_elems("output/elems5.vtk");


    mesher3.mark_faces(mesh3, surf, &simucell);
    mesher3.mark_elems_byvol(mesh3, &simucell);
//    mesh3->recalc("pqAa"+Vmax);
    mesh3->write_faces("output/faces6.vtk");
    mesh3->write_elems("output/elems6.vtk");

    end_msg(t0);

    start_msg("======= Femocs finished! =======\n");
}

} /* namespace femocs */


// ================================================
// ==== Main function to run standalone Femocs ====
// ************************************************
#if not LIBARYMODE

int main(int argc, char* argv[]) {
    using namespace femocs;

    Femocs femocs("/path/to/input.script");

    // If input file specified on command line, use that instead of default
    if (argc >= 2) femocs.conf.infile = argv[1];
    // The same with number of OpenMP threads
    if (argc >= 3) femocs.conf.nt = stod(argv[2]);


    // Make some dummy variables that resemble Helmod input/output format
    const int N = 2;
    int i, j;
    double E0 = 0.0;
    double grid_spacing[3] = {femocs.conf.latconst, femocs.conf.latconst, femocs.conf.latconst};
    double*** BC = new double**[N];
    double*** phi_guess = new double**[N];

    for (i = 0; i < N; ++i) {
        BC[i] = new double*[N];
        phi_guess[i] = new double*[N];
        for (j = 0; j < N; ++j) {
            BC[i][j] = new double[N];
            phi_guess[i][j] = new double[N];
        }
    }

    // Run the actual script
    femocs.run_femocs(E0, BC, phi_guess, grid_spacing);

    return 0;
}

#endif
