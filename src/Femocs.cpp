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
#include "Coarseners.h"
#include "SolutionReader.h"

using namespace std;

// Femocs constructor, specifies simulation parameters
Femocs::Femocs(string message) : solution_valid(false) {
    start_msg(double t0, "======= Femocs started! =======\n");
    conf.external_msg = message;

//    conf.infile = "input/rough111.ckx";
//    conf.infile = "input/mushroom2.ckx";

//    conf.infile = "input/tower_hl2p5.ckx";
//    conf.infile = "input/nanotip_hr5.ckx";
//    conf.latconst = 2.0;

    conf.infile = "input/nanotip_medium.xyz";
    conf.latconst = 3.61;

    conf.coord_cutoff = 3.1;         // coordination analysis cut-off radius

    conf.nnn = 12;                   // number of nearest neighbours in bulk
    conf.mesher = "tetgen";          // mesher algorithm
    conf.mesh_quality = "2.0";
    conf.nt = 4;                     // number of OpenMP threads
    conf.radius = 30.0;              // inner radius of coarsening cylinder
    conf.coarse_factor = 0.5;        // coarsening factor; bigger number gives coarser surface
    conf.smooth_factor = 0.6;        // surface smoothing factor; bigger number gives smoother surface
    conf.n_bins = 20;                // number of bins in histogram smoother
    conf.postprocess_marking = true; // make extra effort to mark correctly the vacuum nodes in shadow area
    conf.rmin_rectancularize = conf.latconst / 1.0; // 1.5+ for <110> simubox, 1.0 for all others
    conf.movavg_width = 1.0;         // width of moving average in electric field smoother

    conf.refine_apex = false;        // refine nanotip apex
    conf.significant_distance = conf.latconst / 2.0;
}

// Destructor, deletes data and prints bye-bye-message
Femocs::~Femocs() {
    start_msg(double t0, "======= Femocs finished! =======\n");
}

// Workhorse function to run Femocs simulation
const void Femocs::run(double E_field, string message) {
    double t0, tstart;  // Variables used to measure the start time of a section
    bool success;
    solution_valid = false;

    conf.neumann = E_field;
    tstart = omp_get_wtime();

//    start_msg(t0, "=== Outputting AtomReader...");
//    reader.output("output/reader.xyz");
//    end_msg(t0);

    start_msg(t0, "=== Comparing with previous run...");
    if ( reader.equals_previous_run(conf.significant_distance) ) return;
    end_msg(t0);

    // ====================================
    // ===== Converting imported data =====
    // ====================================

    start_msg(t0, "=== Extracting surface...");
    femocs::Surface dense_surf(conf.latconst, conf.nnn);
    dense_surf.extract_surface(&reader);
    dense_surf.output("output/surface_dense.xyz");
    end_msg(t0);

    start_msg(t0, "=== Coarsening surface...");
    femocs::Surface coarse_surf(conf.latconst, conf.nnn);
    femocs::Coarseners coarseners;
    coarseners.generate(dense_surf, conf.radius, conf.coarse_factor);

    if (conf.coarse_factor > 0)
        coarse_surf = dense_surf.coarsen(coarseners, &reader.sizes);
    else
        coarse_surf = dense_surf.rectangularize(&reader.sizes, conf.rmin_rectancularize);

    coarse_surf.output("output/surface_coarse.xyz");
    end_msg(t0);

    start_msg(t0, "=== Smoothing surface...");
    coarse_surf.smoothen(conf.radius, conf.smooth_factor, 2 * conf.coord_cutoff);
    coarse_surf.output("output/surface_smooth.xyz");
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
    bulk.generate_simple(&reader.sizes);
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

    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality
    success = mesher.generate_mesh(bulk, coarse_surf, vacuum, "rnQq" + conf.mesh_quality);
    check_success(success, "Triangulation failed! Field calculation will be skipped!");

    big_mesh.write_nodes("output/nodes_generated.xyz");
    big_mesh.write_edges("output/edges_generated.vtk");
    big_mesh.write_faces("output/faces_generated.vtk");
    big_mesh.write_elems("output/elems_generated.vtk");
    end_msg(t0);

    start_msg(t0, "=== Making big mesh appendices...");
    mesher.generate_mesh_appendices();
    big_mesh.write_edges("output/edges_appended.vtk");
    big_mesh.write_faces("output/faces_appended.vtk");
    end_msg(t0);

    start_msg(t0, "=== Marking big mesh...");
    success = mesher.mark_mesh(conf.postprocess_marking);
    big_mesh.write_nodes("output/nodes_marked.xyz");
    big_mesh.write_edges("output/edges_marked.vtk");
    big_mesh.write_faces("output/faces_marked.vtk");
    big_mesh.write_elems("output/elems_marked.vtk");
    check_success(success, "Mesh marking failed! Field calcualtion will be skipped!");
    end_msg(t0);

    start_msg(t0, "=== Separating vacuum and bulk mesh...");
    femocs::Mesh vacuum_mesh(conf.mesher);
    femocs::Mesh bulk_mesh(conf.mesher);

    mesher.separate_meshes_vol2(&vacuum_mesh, &bulk_mesh, "rnQ");

    bulk_mesh.write_nodes("output/nodes_bulk.xyz");
    bulk_mesh.write_edges("output/edges_bulk.vtk");
    bulk_mesh.write_faces("output/faces_bulk.vtk");
    bulk_mesh.write_elems("output/elems_bulk.vtk");

    vacuum_mesh.write_nodes("output/nodes_vacuum.xyz");
    vacuum_mesh.write_edges("output/edges_vacuum.vtk");
    vacuum_mesh.write_faces("output/faces_vacuum.vtk");
    vacuum_mesh.write_elems("output/elems_vacuum.vtk");
    end_msg(t0);


    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    tethex::Mesh tethex_big;
    tethex_big.read_femocs(&big_mesh, {TYPES.NONE, TYPES.BULK, TYPES.SURFACE, TYPES.VACUUM});
    tethex_big.convert();
    end_msg(t0);

    start_msg(t0, "=== Smoothing hexahedra...");
    tethex_big.smoothen(conf.radius, conf.smooth_factor, conf.coord_cutoff);
    tethex_big.export_vertices(&vacuum_mesh);
    end_msg(t0);

    start_msg(t0, "=== Separating hexahedra...");
    tethex::Mesh tethex_bulk;
    tethex::Mesh tethex_vacuum;
    tethex_big.separate_meshes(&tethex_bulk, &tethex_vacuum);

    tethex_bulk.calc_hex_qualities();
    tethex_vacuum.calc_hex_qualities();

    end_msg(t0);

    // ==============================
    // ===== Running FEM solver =====
    // ==============================

    femocs::DealII laplace;
    laplace.set_neumann(conf.neumann);

    start_msg(t0, "=== Importing tethex mesh into Deal.II...");
    success = laplace.import_mesh_wo_faces(&tethex_vacuum);
    check_success(success, "Importing mesh to Deal.II failed! Field calculation will be skipped!");
    end_msg(t0);

    if (conf.refine_apex) {
        start_msg(t0, "=== Refining mesh...");
        bulk.calc_statistics();
        // TODO Fix invalid bulk zmax
        femocs::Point3 origin((bulk.sizes.xmin+bulk.sizes.xmax)/2, (bulk.sizes.ymin+bulk.sizes.ymax)/2, 34);
        laplace.smooth_and_refine_mesh(origin, 7*conf.latconst);
        laplace.output_mesh("output/elems_dealii.vtk");
        end_msg(t0);
    }

    start_msg(t0, "=== System setup...");
    laplace.setup_system();
    end_msg(t0);

#if VERBOSE
    cout << "Vacuum: #elems=" << vacuum_mesh.get_n_elems() << ",\t#faces="
           << vacuum_mesh.get_n_faces() << ",\t#nodes=" << vacuum_mesh.get_n_nodes() << endl;
    cout << "Bulk:   #elems=" << bulk_mesh.get_n_elems() << ",\t#faces=" << bulk_mesh.get_n_faces()
           << ",\t#nodes=" << bulk_mesh.get_n_nodes() << endl;
    cout << "Tethex: #elems=" << tethex_vacuum.get_n_hexahedra() << ",\t#faces="
           << tethex_vacuum.get_n_quadrangles() << ",\t#nodes=" << tethex_vacuum.get_n_vertices() << endl;
    cout << "# degrees of freedom: " << laplace.get_n_dofs() << endl;
    cout << "# input atoms: " << reader.get_n_atoms() << endl;
#endif

    start_msg(t0, "=== Marking the boundary...");
    laplace.mark_boundary_faces(&reader.sizes);
    end_msg(t0);

    start_msg(t0, "=== Assembling system...");
    laplace.assemble_system();
    end_msg(t0);

    start_msg(t0, "=== Solving...");
    laplace.solve_cg();
//    laplace.solve_umfpack();
    end_msg(t0);

    // ========================================
    // ===== Post-processing FEM solution =====
    // ========================================

    start_msg(t0, "=== Extracting solution...");
    femocs::Medium medium = vacuum_mesh.to_medium();
    femocs::SolutionReader solution(&vacuum_mesh);
    solution.extract_solution(laplace);
    end_msg(t0);

    start_msg(t0, "=== Cleaning solution...");
    solution.clean(conf.n_bins);
    end_msg(t0);
    solution.print_statistics();

    start_msg(t0, "=== Interpolating solution...");
    dense_surf.sort_atoms(0, 1, "up");
    interpolation.extract_interpolation(&vacuum_mesh, &solution, dense_surf);
    end_msg(t0);

    start_msg(t0, "=== Cleaning interpolation...");
    interpolation.clean(conf.n_bins);
    end_msg(t0);

    start_msg(t0, "=== Saving results...");
    solution.output("output/results.xyz");
    interpolation.output("output/interpolation" + message + ".xyz");
    laplace.output_results("output/results.vtk");

    coarseners.write("output/coarseners" + message + ".vtk");

    tethex_bulk.write_vtk_elems("output/bulk_smooth" + message + ".vtk");
    tethex_vacuum.write_vtk_elems("output/vacuum_smooth" + message + ".vtk");
    end_msg(t0);

    if (conf.significant_distance > 0.0) reader.save_current_run_points();

    cout << "\nTotal time of Femocs: " << omp_get_wtime() - tstart << "\n";
    solution_valid = true;
}

const void Femocs::import_atoms(string file_name) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    if (file_name == "")
        reader.import_file(conf.infile);
    else
        reader.import_file(file_name);
    end_msg(t0);

    start_msg(t0, "=== Finding coordinations and atom types...");
    reader.calc_coordination(conf.nnn, conf.coord_cutoff);
    reader.extract_types(conf.nnn, conf.latconst);
    end_msg(t0);
}

const void Femocs::import_atoms(int n_atoms, const double* coordinates, const double* box, const int* nborlist) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    reader.import_parcas(n_atoms, coordinates, box);
    end_msg(t0);

    start_msg(t0, "=== Finding coordinations and atom types...");
    reader.calc_coordination(conf.nnn, conf.coord_cutoff, nborlist);
    reader.extract_types(conf.nnn, conf.latconst);
    end_msg(t0);
}

const void Femocs::import_atoms(int n_atoms, double* x, double* y, double* z, int* types) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);

    start_msg(t0, "=== Finding coordinations and atom types...");
    reader.calc_coordination(conf.nnn, conf.coord_cutoff);
    reader.extract_types(conf.nnn, conf.latconst);
    end_msg(t0);
}

const void Femocs::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    if (!solution_valid) return;
    start_msg(double t0, "=== Exporting results...");
    interpolation.export_helmod(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);
}

const void Femocs::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm, const int* nborlist) {
    if (!solution_valid) return;
    start_msg(double t0, "=== Exporting results...");
    interpolation.export_helmod(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);
}

const void femocs_speaker(string path) {
}
