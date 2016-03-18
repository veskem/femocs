/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include <omp.h>
#include <iostream>
#include <memory>

#include "Femocs.h"
#include "Mesher.h"
#include "SurfaceExtractor.h"
#include "Vacuum.h"
#include "tethex.h"
#include "DealII.h"
#include "Mesh.h"

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
//*
//    conf.infile = "input/rough110.ckx";
    conf.infile = "input/mushroom1.ckx";
    conf.latconst = 2.0;        // lattice constant
    conf.coord_cutoff = 3.3;    // coordination analysis cut off radius
//*/
/*
    conf.infile = "input/nanotip_medium.xyz";
    conf.latconst = 3.61;       // lattice constant
    conf.coord_cutoff = 3.4;	// coordination analysis cut off radius
//*/

    conf.nnn = 12;						    // number of nearest neighbours in bulk
    conf.extracter = "coordination";	    // surface extraction algorithm
    conf.mesher = "tetgen";				    // mesher algorithm
    conf.nt = 4;						    // number of OpenMP threads
    conf.poly_degree = 1;                   // finite element polynomial degree
    conf.neumann = 10.0;                    // Neumann boundary condition value
    
    return conf;
}

// Initialise general data about simulation cell
const Femocs::SimuCell Femocs::init_simucell() const {
    SimuCell sc;
    sc.xmin = sc.xmax = sc.ymin = sc.ymax = sc.zmin = sc.zmax = 0;
    sc.zminbox = 0.0;
    sc.zmaxbox = conf.latconst*20; // Electric field will be applied that much higer from the highest coordinate of simubox
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
    double t0, tstart;
    // Max tetrahedron volume is ~100000x the volume of regular tetrahedron with edge == latconst
    //string Vmax = to_string(  (int)(100000.0*0.118*pow(conf.latconst,3.0) ) );

    tstart = start_msg("======= Femocs started =======");
        
#if not LIBRARYMODE
    t0 = start_msg("=== Reading atoms from " + conf.infile);
    reader.import_file(conf.infile, &simucell);
    end_msg(t0);

    t0 = start_msg("=== Calculating coordinations of atoms...");
    reader.calc_coordination(&simucell, conf.coord_cutoff, conf.nnn);
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
    simucell.zmin = surf->sizes.zmin - surf->sizes.latconst;
    simucell.zbox = simucell.zmax - simucell.zmin;
    end_msg(t0);

    t0 = start_msg("=== Extracting bulk...");
    SurfaceExtractor bulk_extractor(&conf);
//    auto bulk = bulk_extractor.extract_bulk(&reader.data, &simucell);
    auto bulk = bulk_extractor.extract_truncated_bulk(&reader.data, &simucell);
    bulk_extractor.rectangularize(bulk, &simucell);
    end_msg(t0);

    t0 = start_msg("=== Generating vacuum...");
    Vacuum vacuum;
    vacuum.generate_simple(&simucell);
    end_msg(t0);

    t0 = start_msg("=== Writing surface, bulk and vacuum to output...");
    surf->output("output/surface.xyz");
    bulk->output("output/bulk.xyz");
    vacuum.output("output/vacuum.xyz");
    end_msg(t0);

// ===============================================
    Mesher mesher(conf.mesher, bulk->sizes.latconst);
    shared_ptr<Mesh> big_mesh(new Mesh());
    shared_ptr<Mesh> bulk_mesh(new Mesh());
    shared_ptr<Mesh> vacuum_mesh(new Mesh());
    
    t0 = start_msg("=== Making big mesh...");
    mesher.get_volume_mesh(big_mesh, bulk, &vacuum, "Q");
    big_mesh->write_faces("output/faces0.vtk");
    big_mesh->write_elems("output/elems0.vtk");

    big_mesh->recalc("rq2.514"); //a"+Vmax);
    big_mesh->write_faces("output/faces1.vtk");
    big_mesh->write_elems("output/elems1.vtk");
    end_msg(t0);

    t0 = start_msg("=== Separating vacuum and bulk mesh...");
//    mesher.extract_vacuum_mesh(vacuum_mesh, big_mesh, bulk->get_n_atoms(), surf->get_n_atoms(), simucell.zmin, "rQ");
//    mesher.extract_bulk_mesh(bulk_mesh, big_mesh, bulk->get_n_atoms(), surf->get_n_atoms(), simucell.zmin, "rQ");
    mesher.separate_meshes(bulk_mesh, vacuum_mesh, big_mesh, bulk->get_n_atoms(), surf->get_n_atoms(), simucell.zmin, "rQ");

    bulk_mesh->write_faces("output/faces_bulk.vtk");
    bulk_mesh->write_elems("output/elems_bulk.vtk");
    vacuum_mesh->write_faces("output/faces_vacuum.vtk");
    vacuum_mesh->write_elems("output/elems_vacuum.vtk");
    end_msg(t0);

//    t0 = start_msg("Making test mesh...");
//    Mesher simplemesher(conf.mesher, bulk->sizes.latconst);
//    auto testmesh = simplemesher.get_simple_mesh();
//    testmesh->recalc("rQ");
//    testmesh->write_elems("output/testelems.vtk");
//    testmesh->write_faces("output/testfaces.vtk");
//    end_msg(t0);

    t0 = start_msg("Converting tetrahedra to hexahedra...");
    tethex::Mesh tethex_mesh;
    tethex_mesh.read_femocs(vacuum_mesh);
    tethex_mesh.convert();
    end_msg(t0);

//    t0 = start_msg("Writing tethex to file...");
//    tethex_mesh.write("output/tethex.msh");
//    tethex_mesh.write_vtk_faces("output/tethex_faces.vtk");
//    tethex_mesh.write_vtk_elems("output/tethex_elems.vtk");
//    end_msg(t0);
    
    DealII laplace(conf.poly_degree, conf.neumann, &simucell);
        
    t0 = start_msg("Importing tethex mesh into Deal.II...");
	laplace.import_tethex_mesh(&tethex_mesh);
    end_msg(t0);
    
    t0 = start_msg("System setup...");
    laplace.setup_system();
    end_msg(t0);

    t0 = start_msg("Marking the boundary...");
    laplace.mark_boundary();
    end_msg(t0);
    
    t0 = start_msg("Assembling system...");
    laplace.assemble_system();
    end_msg(t0);

    t0 = start_msg("Solving...");
    laplace.solve_cg();
//    laplace.solve_umfpack();
    end_msg(t0);

    t0 = start_msg("Outputting results...");
    laplace.output_results("output/final-results.vtk");
    end_msg(t0);

    cout << "\nTotal time: " << omp_get_wtime() - tstart << "\n";
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
