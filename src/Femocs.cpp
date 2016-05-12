/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"

#include "DealII.h"
#include "Media.h"
#include "Mesh.h"
#include "Mesher.h"
#include "Tethex.h"

#include <omp.h>

using namespace std;

// Femocs constructor, specifies simulation parameters
Femocs::Femocs(string file_name) {
    start_msg(double t0, "\n======= Femocs started! =======\n");
    //*
    //conf.infile = "input/rough111.ckx";
    conf.infile = "input/mushroom2.ckx";
    conf.latconst = 2.0;        // lattice constant
    conf.coord_cutoff = 0.0; //3.3;    // coordination analysis cut off radius
    //*/
    /*
     conf.infile = "input/nanotip_big.xyz";
     conf.latconst = 3.61;       // lattice constant
     conf.coord_cutoff = 0.0; //3.4;   // coordination analysis cut off radius
     //*/

    conf.nnn = 12;                          // number of nearest neighbours in bulk
    conf.mesher = "tetgen";                 // mesher algorithm
    conf.nt = 4;                            // number of OpenMP threads
    conf.poly_degree = 1;                   // finite element polynomial degree
}

// Destructor, deletes data and prints bye-bye-message
Femocs::~Femocs() {
    start_msg(double t0, "\n======= Femocs finished! =======\n");
}

// Workhorse function to run Femocs simulation
const void Femocs::run(double E_field, double*** phi) {

    double t0, tstart;  // Variables used to measure the start time of a section
    tstart = omp_get_wtime();
    conf.neumann = E_field;

    // ===== Converting imported data =====

    start_msg(t0, "=== Extracting surface...");
    femocs::Surface surf(conf.latconst, conf.nnn);
    surf.extract_surface(&reader);
    surf.rectangularize(&reader.sizes);
    surf.output("output/surface.xyz");
    end_msg(t0);

    start_msg(t0, "=== Resizing simulation box...");
    // Electric field is applied 20 lattice constants above the surface highest point
    if (surf.get_n_atoms() > 0)
        reader.resize_box(surf.sizes.zmin - conf.latconst, surf.sizes.zmax + 20 * conf.latconst);
    else
        reader.resize_box(reader.sizes.zmin, reader.sizes.zmax + 20 * conf.latconst);
    end_msg(t0);

    start_msg(t0, "=== Extracting bulk...");
    femocs::Bulk bulk(conf.latconst, conf.nnn);
//    bulk.extract_bulk(&reader);
//    bulk.extract_truncated_bulk(&reader);
//    bulk.rectangularize(&reader.sizes);
    bulk.generate_simple(&reader.sizes, &reader.types);
    bulk.output("output/bulk.xyz");
    end_msg(t0);

    start_msg(t0, "=== Generating vacuum...");
    femocs::Vacuum vacuum;
    vacuum.generate_simple(&reader.sizes);
    vacuum.output("output/vacuum.xyz");
    end_msg(t0);

    // ===== Making FEM mesh =====

    femocs::Mesher mesher(conf.mesher, conf.latconst);

    start_msg(t0, "=== Making big mesh...");
    femocs::Mesh big_mesh(conf.mesher);

    mesher.get_volume_mesh(&big_mesh, &bulk, &surf, &vacuum, "Q");
    big_mesh.write_faces("output/faces_0.vtk");
    big_mesh.write_elems("output/elems_0.vtk");
    big_mesh.recalc("rq2.514");
    big_mesh.write_faces("output/faces_1.vtk");
    big_mesh.write_elems("output/elems_1.vtk");
    //mesher.generate_monolayer_surf_faces(&big_mesh, bulk.get_n_atoms(), surf.get_n_atoms());
    mesher.generate_surf_faces(&big_mesh, surf.get_n_atoms());
    big_mesh.write_faces("output/faces_2.vtk");
    big_mesh.write_elems("output/elems_2.vtk");
    end_msg(t0);

    start_msg(t0, "=== Marking nodes...");
    mesher.mark_nodes(&big_mesh, &reader.types, surf.get_n_atoms());
    big_mesh.write_nodes("output/nodes.xyz");
    end_msg(t0);

    start_msg(t0, "=== Separating vacuum and bulk mesh...");
    femocs::Mesh bulk_mesh(conf.mesher);
    femocs::Mesh vacuum_mesh(conf.mesher);

//    mesher.separate_vacuum_mesh(&vacuum_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), bulk.sizes.zmin, "rQ");
//    mesher.separate_bulk_mesh(&bulk_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), bulk.sizes.zmin, "rQ");
//    mesher.separate_meshes(&bulk_mesh, &vacuum_mesh, &big_mesh, bulk.get_n_atoms(), surf.get_n_atoms(), reader.sizes.zminbox, "rQ");

    mesher.separate_meshes_bymarker(&bulk_mesh, &vacuum_mesh, &big_mesh, &reader.types, "rQ");

    bulk_mesh.write_faces("output/faces_bulk.vtk");
    bulk_mesh.write_elems("output/elems_bulk.vtk");
    vacuum_mesh.write_faces("output/faces_vacuum.vtk");
    vacuum_mesh.write_elems("output/elems_vacuum.vtk");
    vacuum_mesh.write_nodes("output/nodes_vacuum.xyz");
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

#if DEBUGMODE
    cout << "Bulk:   #elems=" << bulk_mesh.get_n_elems() << ", #faces=" << bulk_mesh.get_n_faces() << endl;
    cout << "Vacuum: #elems=" << vacuum_mesh.get_n_elems() << ", #faces=" << vacuum_mesh.get_n_faces() << endl;
    cout << "Tetgen: #elems=" << tethex_mesh.get_n_hexahedra() << ", #faces=" << tethex_mesh.get_n_quadrangles() << endl;
#endif

//    start_msg(t0, "Writing tethex to file...");
////    tethex_mesh.write("output/tethex.msh");
//    tethex_mesh.write_vtk_faces("output/tethex_faces.vtk");
//    tethex_mesh.write_vtk_elems("output/tethex_elems.vtk");
//    end_msg(t0);

    // ===== Running FEM solver =====

    femocs::DealII laplace(conf.poly_degree, conf.neumann);

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

    start_msg(t0, "=== Extracting solution...");
    laplace.extract_solution_at_medium(surf);
    end_msg(t0);

    start_msg(t0, "=== Outputting results...");
    laplace.output_results("output/results_vacuum.vtk");
    laplace.output_results("output/results_surface.xyz");
    end_msg(t0);

    cout << "\nTotal time: " << omp_get_wtime() - tstart << "\n";
}

const void Femocs::import_atoms(int n_atoms, double* x, double* y, double* z, int* types) {
    double t0;

    start_msg(t0, "=== Importing atoms...");

    if (n_atoms < 1)
        reader.import_file(conf.infile);
    else
        reader.import_helmod(n_atoms, x, y, z, types);

    end_msg(t0);

    start_msg(t0, "=== Outputting AtomReader atoms...");
    reader.output("output/atoms.xyz");
    end_msg(t0);

    start_msg(t0, "=== Calculating coordinations of atoms...");
    reader.calc_coordination(conf.coord_cutoff, conf.nnn);
    end_msg(t0);


}

const void femocs_speaker(string path) {
    double x[5] = { 1, 2, 3, 4, 5 };
    Femocs f(path);
}
