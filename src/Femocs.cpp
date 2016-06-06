/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include <omp.h>
#include "Femocs.h"

#include "Macros.h"
#include "Medium.h"
#include "Media.h"
#include "Mesh.h"
#include "Mesher.h"
#include "Tethex.h"
#include "DealII.h"

using namespace std;

// Femocs constructor, specifies simulation parameters
Femocs::Femocs(string file_name) {
    start_msg(double t0, "======= Femocs started! =======\n");
    //*
    //conf.infile = "input/rough111.ckx";
    conf.infile = "input/mushroom1.ckx";
    conf.latconst = 2.0;        // lattice constant
    conf.coord_cutoff = 0.0; //3.3;    // coordination analysis cut off radius
    //*/
    /*
     conf.infile = "input/nanotip_big.xyz";
     conf.latconst = 3.61;       // lattice constant
     conf.coord_cutoff = 0.0; //3.4;   // coordination analysis cut off radius
     //*/

    conf.nnn = 12;                  // number of nearest neighbours in bulk
    conf.mesher = "tetgen";         // mesher algorithm
    conf.mesh_quality = "2.414";
    conf.nt = 4;                    // number of OpenMP threads
    conf.rmin_coarse = 16.5;        // inner radius of coarsening cylinder
    conf.coarse_factor = 0.5;       // coarsening factor; bigger number gives coarser surface
    conf.postprocess_marking = true; // make extra effort to mark correctly the vacuum nodes in shadow area
    conf.rmin_rectancularize = conf.latconst / 1.0; // 1.5+ for <110> simubox, 1.0 for all others
}

// Destructor, deletes data and prints bye-bye-message
Femocs::~Femocs() {
    start_msg(double t0, "======= Femocs finished! =======\n");
}

// Workhorse function to run Femocs simulation
const void Femocs::run(double E_field) {
    double t0, tstart;  // Variables used to measure the start time of a section
    tstart = omp_get_wtime();

    conf.neumann = E_field;

    // ====================================
    // ===== Converting imported data =====
    // ====================================

    femocs::Surface dense_surf(conf.latconst, conf.nnn);
    femocs::Surface coarse_surf(conf.latconst, conf.nnn);

    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract_surface(&reader);

    if (conf.coarse_factor > 0)
        coarse_surf = dense_surf.coarsen(conf.rmin_coarse, conf.coarse_factor, &reader.sizes);
    else
        coarse_surf = dense_surf.rectangularize(&reader.sizes, conf.rmin_rectancularize);

    dense_surf.output("output/dense_surface.xyz");
    coarse_surf.output("output/coarse_surface.xyz");
    end_msg(t0);

    start_msg(t0, "=== Resizing simulation box...");
    // Electric field is applied 100 lattice constants above the highest point of surface
    // and bulk is extended 20 lattice constants below the minimum point of surface
    const double zmaxbox = coarse_surf.sizes.zmax + 100 * conf.latconst;
    const double zminbox = coarse_surf.sizes.zmin - 20 * conf.latconst;
    reader.resize_box(zminbox, zmaxbox);
    end_msg(t0);

    start_msg(t0, "=== Extracting bulk...");
    femocs::Bulk bulk(conf.latconst, conf.nnn);
    //    bulk.extract_bulk(&reader);
    //    bulk.extract_truncated_bulk(&reader);
    bulk.generate_simple(&reader.sizes, &reader.types);
    bulk.output("output/bulk.xyz");
    end_msg(t0);

    start_msg(t0, "=== Generating vacuum...");
    femocs::Vacuum vacuum;
    vacuum.generate_simple(&reader.sizes);
    vacuum.output("output/vacuum.xyz");
    end_msg(t0);

    // ===========================
    // ===== Making FEM mesh =====
    // ===========================

    start_msg(t0, "=== Making big mesh...");
    femocs::Mesh big_mesh(conf.mesher);
    femocs::Mesher mesher(&big_mesh);

    mesher.generate_mesh(bulk, coarse_surf, vacuum, "rQq" + conf.mesh_quality);
    big_mesh.write_nodes("output/nodes_generated.xyz");
    big_mesh.write_faces("output/faces_generated.vtk");
    big_mesh.write_elems("output/elems_generated.vtk");
    end_msg(t0);

    start_msg(t0, "=== Marking nodes...");
    mesher.mark_mesh(&reader.types, conf.postprocess_marking);
    big_mesh.write_nodes("output/nodes_big.xyz");
    big_mesh.write_faces("output/faces_big.vtk");
    big_mesh.write_elems("output/elems_big.vtk");
    end_msg(t0);

    start_msg(t0, "=== Separating vacuum and bulk mesh...");
    femocs::Mesh bulk_mesh(conf.mesher);
    femocs::Mesh vacuum_mesh(conf.mesher);

    mesher.separate_meshes(&bulk_mesh, &vacuum_mesh, &reader.types, "rQ");
    bulk_mesh.write_nodes("output/nodes_bulk.xyz");
    bulk_mesh.write_faces("output/faces_bulk.vtk");
    bulk_mesh.write_elems("output/elems_bulk.vtk");
    vacuum_mesh.write_nodes("output/nodes_vacuum.xyz");
    vacuum_mesh.write_faces("output/faces_vacuum.vtk");
    vacuum_mesh.write_elems("output/elems_vacuum.vtk");
    end_msg(t0);

//     start_msg(t0, "=== Making test mesh...");
//     femocs::Mesh testmesh(conf.mesher);
//     testmesh.generate_simple();
//     testmesh.write_elems("output/testelems.vtk");
//     testmesh.write_faces("output/testfaces.vtk");
//     end_msg(t0);

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    tethex::Mesh tethex_mesh;
    tethex_mesh.read_femocs(&vacuum_mesh);
    tethex_mesh.convert();
    end_msg(t0);

//     start_msg(t0, "=== Writing tethex to file...");
//     tethex_mesh.write("output/tethex.msh");
//     tethex_mesh.write_vtk_faces("output/tethex_faces.vtk");
//     tethex_mesh.write_vtk_elems("output/tethex_elems.vtk");
//     end_msg(t0);

    // ==============================
    // ===== Running FEM solver =====
    // ==============================
    femocs::DealII laplace;
    laplace.set_neumann(conf.neumann);

    start_msg(t0, "=== Importing tethex mesh into Deal.II...");
    laplace.import_tethex_mesh(&tethex_mesh);
    end_msg(t0);

    start_msg(t0, "=== System setup...");
    laplace.setup_system();
    end_msg(t0);

#if DEBUGMODE
    cout << "Bulk:   #elems=" << bulk_mesh.get_n_elems() << ", #faces=" << bulk_mesh.get_n_faces() << endl;
    cout << "Vacuum: #elems=" << vacuum_mesh.get_n_elems() << ", #faces=" << vacuum_mesh.get_n_faces() << endl;
    cout << "Tetgen: #elems=" << tethex_mesh.get_n_hexahedra() << ", #faces=" << tethex_mesh.get_n_quadrangles() << endl;
    cout << "# degrees of freedom: " << laplace.get_n_dofs() << endl;
    cout << "# input atoms: " << reader.get_n_atoms() << endl;
#endif

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
//     femocs::Medium medium = vacuum_mesh.to_medium();
//     solution.extract_solution_at_medium(&laplace, medium);
    solution.extract_solution(&laplace, coarse_surf);
    end_msg(t0);

    start_msg(t0, "=== Outputting results...");
    solution.output("output/results_surface.xyz");
    end_msg(t0);

#if DEBUGMODE
    cout << "\nTotal time: " << omp_get_wtime() - tstart << "\n";
#endif
}

const void Femocs::import_atoms(int n_atoms, double* x, double* y, double* z, int* types) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    if (n_atoms < 1)
        reader.import_file(conf.infile);
    else
        reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);

    start_msg(t0, "=== Calculating coordinations of atoms...");
    reader.calc_coordination(conf.coord_cutoff, conf.nnn);
    end_msg(t0);
}

const void Femocs::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    start_msg(double t0, "=== Exporting results...");
    solution.export_helmod(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);
}

const void femocs_speaker(string path) {
}
