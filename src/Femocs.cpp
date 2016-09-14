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
Femocs::Femocs(string file_name) :
        solution_valid(false) {
    start_msg(double t0, "======= Femocs started! =======\n");
    conf.path_to_script = file_name;
    //*
//    conf.infile = "input/rough111.ckx";
//    conf.infile = "input/mushroom2.ckx";
    conf.infile = "input/nanotip_hr5.ckx";
    conf.latconst = 2.0;        // lattice constant
    conf.coord_cutoff = 3.1;    // coordination analysis cut off radius
    //*/
    /*
     conf.infile = "input/nanotip_big.xyz";
     conf.latconst = 3.61;       // lattice constant
     conf.coord_cutoff = 0.0; //3.4;   // coordination analysis cut off radius
     //*/

    conf.nnn = 12;                  // number of nearest neighbours in bulk
    conf.mesher = "tetgen";         // mesher algorithm
    conf.mesh_quality = "2.0";//"2.914";
    conf.nt = 4;                    // number of OpenMP threads
    conf.rmin_coarse = 7.0;        // inner radius of coarsening cylinder
    conf.rmax_coarse = 8000.0;        // radius of constant cutoff coarsening cylinder
    conf.coarse_factor = 0.8;       // coarsening factor; bigger number gives coarser surface
    conf.postprocess_marking = true; //true;//false; // make extra effort to mark correctly the vacuum nodes in shadow area
    conf.rmin_rectancularize = conf.latconst / 1.0; // 1.5+ for <110> simubox, 1.0 for all others
    conf.movavg_width = 2;           // width of moving average in electric field smoother
}

// Destructor, deletes data and prints bye-bye-message
Femocs::~Femocs() {
    start_msg(double t0, "======= Femocs finished! =======\n");
}

// Workhorse function to run Femocs simulation
const void Femocs::run(double E_field) {
    double t0, tstart;  // Variables used to measure the start time of a section
    solution_valid = false;

    conf.neumann = E_field;
    tstart = omp_get_wtime();

    // ====================================
    // ===== Converting imported data =====
    // ====================================

    femocs::Surface dense_surf(conf.latconst, conf.nnn);
    femocs::Surface coarse_surf(conf.latconst, conf.nnn);

    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract_surface(&reader);

    if (conf.coarse_factor > 0)
        coarse_surf = dense_surf.coarsen_vol2(conf.coord_cutoff, conf.rmin_coarse, conf.rmax_coarse,
                conf.coarse_factor, &reader.sizes);
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
    try {
        mesher.generate_mesh(bulk, coarse_surf, vacuum, "rQq" + conf.mesh_quality);
    } catch (int e) {
        cout << "!!! Exception " << e << " occured while making mesh! Field calcualtion will be skipped!\n";
        return;
    }
    big_mesh.write_nodes("output/nodes_generated.xyz");
    big_mesh.write_faces("output/faces_generated.vtk");
    big_mesh.write_elems("output/elems_generated.vtk");
    end_msg(t0);

    start_msg(t0, "=== Marking big mesh...");
    mesher.mark_mesh(conf.postprocess_marking, 1.5*coarse_surf.get_roughness());
    big_mesh.write_nodes("output/nodes_marked.xyz");
    big_mesh.write_nodes("output/nodes_marked.vtk");
    big_mesh.write_faces("output/faces_marked.vtk");
    big_mesh.write_elems("output/elems_marked.vtk");
    end_msg(t0);

    start_msg(t0, "=== Separating vacuum and bulk mesh...");
    femocs::Mesh vacuum_mesh(conf.mesher);
    femocs::Mesh bulk_mesh(conf.mesher);

    mesher.separate_meshes_vol2(&vacuum_mesh, &bulk_mesh, "rQ");


//    mesher.separate_meshes_noclean(&vacuum_mesh, "rQ");
//    mesher.separate_meshes_noclean(&vacuum_mesh, &bulk_mesh, "rQ");

    bulk_mesh.write_nodes("output/nodes_bulk.xyz");
    bulk_mesh.write_faces("output/faces_bulk.vtk");
    bulk_mesh.write_elems("output/elems_bulk.vtk");
    vacuum_mesh.write_nodes("output/nodes_vacuum.xyz");
    vacuum_mesh.write_faces("output/faces_vacuum.vtk");
    vacuum_mesh.write_elems(conf.path_to_script);
    end_msg(t0);

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    tethex::Mesh tethex_mesh;
    int material_ID = 10;
    tethex_mesh.read_femocs(&vacuum_mesh, material_ID);
    tethex_mesh.convert();
    end_msg(t0);

/*
    start_msg(t0, "=== Making union mesh...");
    tethex::Mesh tethex_union_mesh;
    femocs::Mesher bulk_mesher(&bulk_mesh);
    femocs::Mesher vacuum_mesher(&vacuum_mesh);
    vector<int> bi = bulk_mesher.get_new_indxs();
    vector<int> vi = vacuum_mesher.get_new_indxs();
        
    tethex_union_mesh.read_femocs(&bulk_mesh, &vacuum_mesh, bi, vi);
    tethex_union_mesh.convert();
    end_msg(t0);

    start_msg(t0, "=== Writing tethen_averagex to file...");
    tethex_union_mesh.write("output/mushroom.msh");
    end_msg(t0);
//*/     

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

#if VERBOSE
    cout << "Vacuum: #elems=" << vacuum_mesh.get_n_elems() << ",\t#faces="
            << vacuum_mesh.get_n_faces() << ",\t#nodes=" << vacuum_mesh.get_n_nodes() << endl;
//    cout << "Bulk:   #elems=" << bulk_mesh.get_n_elems() << ",\t#faces=" << bulk_mesh.get_n_faces()
//             << ",\t#nodes=" << bulk_mesh.get_n_nodes() << endl;
    cout << "Tethex: #elems=" << tethex_mesh.get_n_hexahedra() << ",\t#faces="
            << tethex_mesh.get_n_quadrangles() << ",\t#nodes=" << tethex_mesh.get_n_vertices() << endl;
    cout << "# degrees of freedom: " << laplace.get_n_dofs() << endl;
    cout << "# input atoms: " << reader.get_n_atoms() << endl;
#endif

    start_msg(t0, "=== Marking the boundary...");
    laplace.mark_boundary(&reader.sizes);
    end_msg(t0);

    start_msg(t0, "=== Assembling system...");
    laplace.assemble_system();
    end_msg(t0);

    start_msg(t0, "=== Solving...");
    laplace.solve_cg();
//    laplace.solve_umfpack();
    end_msg(t0);

    start_msg(t0, "=== Extracting solution...");
//    femocs::Medium medium = vacuum_mesh.to_medium();
//    solution.extract_solution(&laplace, medium);
    solution.extract_solution(&laplace, coarse_surf);
    end_msg(t0);

    start_msg(t0, "=== Smoothing solution...");
    solution.smoothen_result(conf.movavg_width);
    end_msg(t0);
    solution.print_statistics();

//    start_msg(t0, "=== Extracting statistics...");
//    solution.extract_statistics(vacuum_mesh);
//    end_msg(t0);

    start_msg(t0, "=== Outputting results...");
    solution.output("output/results.xyz");
//    laplace.output_results("output/results.vtk");
    end_msg(t0);

//#if VERBOSE
    cout << "\nTotal time of Femocs: " << omp_get_wtime() - tstart << "\n";
//#endif
    solution_valid = true;
}

const void Femocs::import_atoms(int n_atoms, const double* coordinates, const double* box, const int* nborlist) {
    double t0;
    start_msg(t0, "=== Importing atoms...");
    reader.import_parcas(n_atoms, coordinates, box);
    end_msg(t0);

    start_msg(t0, "=== Calculating coordinations...");
    reader.calc_coordination(conf.nnn, conf.coord_cutoff, nborlist);
    end_msg(t0);

    start_msg(t0, "=== Extracting atom types...");
    reader.extract_types(conf.nnn);
    end_msg(t0);

//    start_msg(t0, "=== Outputting AtomReader...");
//    reader.output("output/reader.xyz");
//    end_msg(t0);
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
    reader.calc_coordination(conf.nnn);
    end_msg(t0);

//    start_msg(t0, "=== Outputting AtomReader...");
//    reader.output("output/reader.xyz");
//    end_msg(t0);
}

const void Femocs::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    if (solution_valid) {
        start_msg(double t0, "=== Exporting results...");
        solution.export_helmod(n_atoms, Ex, Ey, Ez, Enorm);
        end_msg(t0);
    }
}

const void Femocs::export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm, const int* nborlist) {
    if (solution_valid) {
        start_msg(double t0, "=== Exporting results...");
        solution.export_helmod(n_atoms, Ex, Ey, Ez, Enorm);
        end_msg(t0);
    }
}

const void femocs_speaker(string path) {
}
