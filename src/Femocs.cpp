/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"

#include "AtomReader.h"
#include "DealII.h"
#include "Macros.h"
#include "Media.h"
#include "Mesh.h"
#include "Mesher.h"
#include "Tethex.h"

using namespace std;
namespace femocs {

// Femocs constructor
Femocs::Femocs(const string file_name) {
    conf = parse_input_script(file_name);
}

// Get simulation parameters
const Femocs::Config Femocs::parse_input_script(const string fileName) const {
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
    conf.mesher = "tetgen";				    // mesher algorithm
    conf.nt = 4;						    // number of OpenMP threads
    conf.poly_degree = 1;                   // finite element polynomial degree
    conf.neumann = 10.0;                    // Neumann boundary condition value

    return conf;
}

// Workhorse function to run Femocs simulation
const void Femocs::run_femocs(const double E0, double*** BC, double*** phi_guess,
        const double* grid_spacing) {
    AtomReader reader;

    double t0, tstart;  // Variables used to measure the start time of a section

    start_msg(tstart, "======= Femocs started =======");

    start_msg(t0, "=== Importing atoms...");
#if STANDALONEMODE
    reader.import_file(conf.infile);
#elif HELMODMODE
    reader.import_helmod(E0, BC, phi_guess,...);
#elif KIMOCSMODE
    reader.import_kimocs();
#else
    cout << "Incorrent MODE definitions in Macros.h !" << endl;
    exit(EXIT_FAILURE);
#endif
    end_msg(t0);

    start_msg(t0, "=== Calculating coordinations of atoms...");
    reader.calc_coordination(conf.coord_cutoff, conf.nnn);
    end_msg(t0);

    start_msg(t0, "=== Extracting surface...");
    Surface surf(conf.latconst, conf.nnn);
    surf.extract_surface(&reader);
    end_msg(t0);

    start_msg(t0, "=== Resizing simulation box...");
    // Electric field is applied 20 lattice constants above the surface highest point
    reader.resize_box(surf.sizes.zmin - conf.latconst, surf.sizes.zmax + 20 * conf.latconst);
    end_msg(t0);

    start_msg(t0, "=== Extracting bulk...");
    Bulk bulk(conf.latconst, conf.nnn);
//    bulk.extract_bulk(&reader);
    bulk.extract_truncated_bulk(&reader);
    bulk.rectangularize(&reader.sizes);
    end_msg(t0);

    start_msg(t0, "=== Generating vacuum...");
    Vacuum vacuum;
    vacuum.generate_simple(&reader.sizes);
    end_msg(t0);

    start_msg(t0, "=== Writing surface, bulk and vacuum to output...");
    surf.output("output/surface.xyz");
    bulk.output("output/bulk.xyz");
    vacuum.output("output/vacuum.xyz");
    end_msg(t0);

// ===============================================
    Mesher mesher(conf.mesher, conf.latconst);

    start_msg(t0, "=== Making big mesh...");
    Mesh big_mesh;
    mesher.get_volume_mesh(&big_mesh, &bulk, &vacuum, "Q");

    big_mesh.write_faces("output/faces0.vtk");
    big_mesh.write_elems("output/elems0.vtk");
    big_mesh.recalc("rq2.514");
    big_mesh.write_faces("output/faces1.vtk");
    big_mesh.write_elems("output/elems1.vtk");
    end_msg(t0);

    start_msg(t0, "=== Separating vacuum and bulk mesh...");
    Mesh bulk_mesh;
    Mesh vacuum_mesh;
//    mesher.separate_vacuum_mesh(&vacuum_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), bulk.sizes.zmin, "rQ");
//    mesher.separate_bulk_mesh(&bulk_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), bulk.sizes.zmin, "rQ");
    mesher.separate_meshes(&bulk_mesh, &vacuum_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), reader.sizes.zminbox, "rQ");

    bulk_mesh.write_faces("output/faces_bulk.vtk");
    bulk_mesh.write_elems("output/elems_bulk.vtk");
    vacuum_mesh.write_faces("output/faces_vacuum.vtk");
    vacuum_mesh.write_elems("output/elems_vacuum.vtk");
    end_msg(t0);

//    start_msg(t0, "Making test mesh...");
//    shared_ptr<Mesh> test_mesh(new Mesh());
//    mesher.get_test_mesh(test_mesh);
//    testmesh->recalc("rQ");
//    testmesh->write_elems("output/testelems.vtk");
//    testmesh->write_faces("output/testfaces.vtk");
//    end_msg(t0);

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    tethex::Mesh tethex_mesh;
    tethex_mesh.read_femocs(&vacuum_mesh);
    tethex_mesh.convert();
    end_msg(t0);

//    t0 = start_msg("Writing tethex to file...");
//    tethex_mesh.write("output/tethex.msh");
//    tethex_mesh.write_vtk_faces("output/tethex_faces.vtk");
//    tethex_mesh.write_vtk_elems("output/tethex_elems.vtk");
//    end_msg(t0);

    DealII laplace(conf.poly_degree, conf.neumann);

    start_msg(t0, "=== Importing tethex mesh into Deal.II...");
    laplace.import_tethex_mesh(&tethex_mesh);
    end_msg(t0);

    start_msg(t0, "=== System setup...");
    laplace.setup_system();
    end_msg(t0);

    start_msg(t0, "=== Marking the boundary...");
    laplace.mark_boundary(&reader.sizes, &reader.types);
    end_msg(t0);

    start_msg(t0, "=== Assembling system...");
    laplace.assemble_system(&reader.types);
    end_msg(t0);

    start_msg(t0, "=== Solving...");
    laplace.solve_cg();
//    laplace.solve_umfpack();
    end_msg(t0);

    start_msg(t0, "=== Outputting results...");
    laplace.output_results("output/final-results.vtk");
    end_msg(t0);

    cout << "\nTotal time: " << omp_get_wtime() - tstart << "\n";
    start_msg(t0, "======= Femocs finished! =======\n");
}

} /* namespace femocs */

// ================================================
// ==== Main function to run standalone Femocs ====
// ************************************************
#if STANDALONEMODE

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
    double grid_spacing[3] = { femocs.conf.latconst, femocs.conf.latconst, femocs.conf.latconst };
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

#endif // STANDALONEMODE
