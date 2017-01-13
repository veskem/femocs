/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"
#include "Coarseners.h"
#include "DealII.h"
#include "Macros.h"
#include "Tethex.h"
#include "TetgenMesh.h"

#include <omp.h>

using namespace std;
namespace femocs {

// specify simulation parameters
Femocs::Femocs(string path_to_conf) : skip_calculations(false) {
    double t0;
    start_msg(t0, "======= Femocs started! =======\n");

    start_msg(t0, "=== Reading configuration parameters...");
    conf.read_all(path_to_conf);
    end_msg(t0);
    conf.print_data();

#if HEATINGMODE
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_solver  = &ch_solver1;
    prev_ch_solver = NULL;
#endif

    // Clear the results from previous run
    if (conf.clear_output) system("rm -rf output");
    if (MODES.WRITEFILE) system("mkdir -p output");
}

// delete data and print bye-bye-message
Femocs::~Femocs() {
    start_msg(double t0, "======= Femocs finished! =======\n");
}

const int Femocs::generate_boundary_nodes(Media& bulk, Media& coarse_surf, Media& vacuum) {
    double t0;
    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract(reader, TYPES.SURFACE);
    dense_surf = dense_surf.clean_lonely_atoms(conf.coord_cutoff);
    end_msg(t0);
    dense_surf.write("output/surface_dense.xyz");
        
    Media stretch_surf;
    stretch_surf = dense_surf.stretch(conf.latconst, conf.box_width);
    stretch_surf.write("output/surface_stretch.xyz");

    Coarseners coarseners;
    coarseners.generate(dense_surf, conf.radius, conf.cfactor, conf.latconst);
    coarseners.write("output/coarseners.vtk");
    
    start_msg(t0, "=== Coarsening surface...");
    coarse_surf = stretch_surf.coarsen(coarseners, stretch_surf.sizes);
    end_msg(t0);
    coarse_surf.write("output/surface_nosmooth.xyz");

    start_msg(t0, "=== Smoothing surface...");
    coarse_surf.smoothen(conf.radius, conf.smooth_factor, 3.0*conf.coord_cutoff);
    end_msg(t0);
    coarse_surf.write("output/surface_coarse.xyz");
    
    start_msg(t0, "=== Generating bulk and vacuum...");
    coarse_surf.calc_statistics();  // calculate zmin and zmax for surface
    bulk.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin - conf.bulk_height * conf.latconst);
    vacuum.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin + conf.box_height * coarse_surf.sizes.zbox);
    reader.resize_box(coarse_surf.sizes.xmin, coarse_surf.sizes.xmax, 
        coarse_surf.sizes.ymin, coarse_surf.sizes.ymax,
        bulk.sizes.zmin, vacuum.sizes.zmax);
    end_msg(t0);
    
    bulk.write("output/bulk.xyz");
    vacuum.write("output/vacuum.xyz");
    
    return 0;
}

const int Femocs::generate_meshes(TetgenMesh& tetmesh_bulk, TetgenMesh& tetmesh_vacuum,
    tethex::Mesh& hexmesh_bulk, tethex::Mesh& hexmesh_vacuum) {

    double t0;
    bool fail;
    
    Media bulk, coarse_surf, vacuum;
    fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
    if(fail) return 1;
    
    start_msg(t0, "=== Making big mesh...");
    TetgenMesh tetmesh_big;
    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality
    fail = tetmesh_big.generate(bulk, coarse_surf, vacuum, "rnQq" + conf.mesh_quality);
    check_message(fail, "Triangulation failed! Field calculation will be skipped!");
    end_msg(t0);

    start_msg(t0, "=== Making surface faces...");
    tetmesh_big.generate_appendices();
    end_msg(t0);

    tetmesh_big.faces.write("output/surface_faces.vtk");

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = tetmesh_big.mark_mesh(conf.postprocess_marking);
    tetmesh_big.nodes.write("output/tetmesh_nodes.xyz");
    tetmesh_big.elems.write("output/tetmesh_elems.vtk");
    check_message(fail, "Mesh marking failed! Field calcualtion will be skipped!");
    end_msg(t0);

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    tethex::Mesh hexmesh_big;
    hexmesh_big.read_femocs(tetmesh_big);
    hexmesh_big.convert();
    end_msg(t0);

    start_msg(t0, "=== Smoothing hexahedra...");
    hexmesh_big.smoothen(conf.radius, conf.smooth_factor, 3.0*conf.coord_cutoff);
    hexmesh_big.export_vertices(tetmesh_big);  // correcting the nodes in tetrahedral mesh
    end_msg(t0);

    start_msg(t0, "=== Separating tetrahedral meshes...");
    tetmesh_big.separate_meshes(tetmesh_bulk, tetmesh_vacuum, "rnQ");
    end_msg(t0);

    if (MODES.VERBOSE)
        cout << "Bulk:   " << tetmesh_bulk << "\nVacuum: " << tetmesh_vacuum << endl;

    start_msg(t0, "=== Separating hexahedral meshes...");
    hexmesh_big.separate_meshes(hexmesh_bulk, hexmesh_vacuum);
    end_msg(t0);
    
    return 0;
}

const int Femocs::solve_laplace(TetgenMesh& tetmesh_vacuum, tethex::Mesh& hexmesh_vacuum) {
    double t0;
    bool fail;
    
    DealII laplace;
    laplace.set_neumann(conf.neumann);

    start_msg(t0, "=== Importing mesh into Deal.II...");
    fail = laplace.import_mesh_wo_faces(hexmesh_vacuum);
    check_message(fail, "Importing mesh to Deal.II failed! Field calculation will be skipped!");
    end_msg(t0);

    if (conf.refine_apex) {
        start_msg(t0, "=== Refining mesh in Deal.II...");
        Point3 origin(dense_surf.sizes.xmid, dense_surf.sizes.ymid, dense_surf.sizes.zmax);
        laplace.refine_mesh(origin, 7*conf.latconst);
        laplace.write_mesh("output/hexmesh_refine.vtk");
        end_msg(t0);
    }

    start_msg(t0, "=== System setup...");
    laplace.setup_system(reader.sizes);
    end_msg(t0);

    start_msg(t0, "=== Assembling system...");
    laplace.assemble_system();
    end_msg(t0);

    start_msg(t0, "=== Solving...");
    laplace.solve_cg();
//    laplace.solve_umfpack();
    end_msg(t0);
    laplace.write("output/result_E_phi.vtk");

    start_msg(t0, "=== Extracting solution...");
    interpolator.extract_solution(laplace, tetmesh_vacuum);
    end_msg(t0);
    interpolator.write("output/result_E_phi.xyz");

    return 0;
}

// workhorse function to generate FEM mesh and to solve differential equation(s)
const int Femocs::run(double elfield, string message) {
    double t0, tstart;  // Variables used to measure the code execution time
    bool fail;

    check_message(skip_calculations, "Atoms haven't moved significantly, " +
            to_string(reader.rms_distance).substr(0,5) + " < " + to_string(conf.distance_tol).substr(0,5)
            + "! Field calculation will be skipped!");

    skip_calculations = true;

    conf.neumann = elfield;
    conf.message = message;
    tstart = omp_get_wtime();

    // ===========================
    // ===== Making FEM mesh =====
    // ===========================

    TetgenMesh tetmesh_bulk, tetmesh_vacuum;
    tethex::Mesh hexmesh_bulk, hexmesh_vacuum;
    
    fail = generate_meshes(tetmesh_bulk, tetmesh_vacuum, hexmesh_bulk, hexmesh_vacuum);
    if(fail) return 1;
    
    tetmesh_bulk.elems.write  ("output/tetmesh_bulk.vtk");
    tetmesh_vacuum.elems.write("output/tetmesh_vacuum.vtk");
    hexmesh_bulk.write_vtk_elems  ("output/hexmesh_bulk" + message + ".vtk");
    hexmesh_vacuum.write_vtk_elems("output/hexmesh_vacuum" + message + ".vtk");

#if HEATINGMODE
    // ====================================
    // ===== Running FEM multi-solver =====
    // ====================================
    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fch::Laplace<3> laplace_solver;
    laplace_solver.import_mesh_directly(hexmesh_vacuum.get_nodes(), hexmesh_vacuum.get_elems());
    end_msg(t0);

    start_msg(t0, "=== Initializing Laplace solver...");
    laplace_solver.set_applied_efield(10.0*elfield);
    laplace_solver.setup_system();
    laplace_solver.assemble_system();
    end_msg(t0);

    if (MODES.VERBOSE) cout << laplace_solver << endl;
        
    start_msg(t0, "=== Running Laplace solver...");
    laplace_solver.solve();
    end_msg(t0);

    if (MODES.WRITEFILE) laplace_solver.output_results("output/result_E_phi" + message + ".vtk");
    
    start_msg(t0, "=== Initializing rho & T solver...");
    ch_solver->reinitialize(&laplace_solver, prev_ch_solver);
    end_msg(t0);
    
    start_msg(t0, "=== Importing mesh to rho & T solver...");
    ch_solver->import_mesh_directly(hexmesh_bulk.get_nodes(), hexmesh_bulk.get_elems());
    ch_solver->setup_system();
    end_msg(t0);

    if (MODES.VERBOSE) cout << *(ch_solver) << endl;
    
    start_msg(t0, "=== Running rho & T solver...\n");
    double temp_error = ch_solver->run_specific(conf.t_error, conf.n_newton, false, "", MODES.VERBOSE, 2.0);
    end_msg(t0);

    if (MODES.WRITEFILE) ch_solver->output_results("output/result_rho_T" + message + ".vtk");
    check_message(temp_error > conf.t_error, "Temperature didn't converge, err=" + to_string(temp_error)) + "! Using previous solution!";
       
    start_msg(t0, "=== Extracting vacuum solution...");
    vacuum_interpolator.extract_solution(&laplace_solver, tetmesh_vacuum);
    end_msg(t0);

    vacuum_interpolator.write("output/result_E_phi.xyz");
    
    start_msg(t0, "=== Extracting bulk solution...");
    bulk_interpolator.extract_solution(ch_solver, tetmesh_bulk);
    end_msg(t0);

    bulk_interpolator.write("output/result_rho_T.xyz");
    
    start_msg(t0, "=== Interpolating E and phi...");
    vacuum_interpolation.interpolate(dense_surf, conf.coord_cutoff*conf.smoothen_solution);
    end_msg(t0);
    
    start_msg(t0, "=== Interpolating T and rho...");
    bulk_interpolation.interpolate(reader, 0);
    end_msg(t0);

    vacuum_interpolation.write("output/interpolation_vacuum" + conf.message + ".vtk");
    bulk_interpolation.write("output/interpolation_bulk" + conf.message + ".vtk");
        
    static bool odd_run = true;
    if (odd_run) {
        ch_solver = &ch_solver2;
        prev_ch_solver = &ch_solver1;
    } else {
        ch_solver = &ch_solver1;
        prev_ch_solver = &ch_solver2;
    }
    odd_run = !odd_run;
#endif

    fail = solve_laplace(tetmesh_vacuum, hexmesh_vacuum);
    if(fail) return 1;

    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.distance_tol);
    end_msg(t0);

    cout << "\nTotal time of Femocs.run: " << omp_get_wtime() - tstart << "\n";
    skip_calculations = false;

    return 0;
}

// import atoms from file
const int Femocs::import_atoms(const string& file_name) {
    double t0;
    string file_type, fname;

    if (file_name == "") fname = conf.infile;
    else fname = file_name;

    file_type = get_file_type(fname);
    expect(file_type == "ckx" || file_type == "xyz", "Unknown file type: " + file_type);

    start_msg(t0, "=== Importing atoms...");
    reader.import_file(fname);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.get_n_atoms() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.get_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations)
        if (file_type == "xyz") {
            start_msg(t0, "=== Calculating coords and atom types...");
            reader.calc_coordination(conf.coord_cutoff);
            reader.extract_types(conf.nnn, conf.latconst);
            end_msg(t0);
        } else {
            start_msg(t0, "=== Calculating coords from atom types...");
            reader.calc_coordination(conf.nnn);
            end_msg(t0);
        }

    return 0;
}

// import atoms from PARCAS
const int Femocs::import_atoms(int n_atoms, double* coordinates, double* box, int* nborlist) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    reader.import_parcas(n_atoms, coordinates, box);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.get_n_atoms() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.get_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating coords and atom types...");
        reader.calc_coordination(conf.nnn, conf.coord_cutoff, nborlist);
        reader.extract_types(conf.nnn, conf.latconst);
        end_msg(t0);
    }

    return 0;
}

// import coordinates and types of atoms
const int Femocs::import_atoms(int n_atoms, double* x, double* y, double* z, int* types) {
    double t0;

    start_msg(t0, "=== Importing atoms...");
    reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.get_n_atoms() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.get_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating coords from atom types...");
        reader.calc_coordination(conf.nnn);
        end_msg(t0);
    }

    return 0;
}

// export the calculated electric field on imported atom coordinates
const int Femocs::export_elfield(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    check_message(interpolator.get_n_atoms() == 0, "No solution to export!");
    double t0;

    if (!skip_calculations) {
        start_msg(t0, "=== Interpolating solution...");
        interpolation.interpolate(dense_surf, conf.smoothen_solution*conf.coord_cutoff);
        end_msg(t0);

        interpolation.write("output/interpolation.movie");
    }

    start_msg(t0, "=== Exporting solution...");
    interpolation.export_solution(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);

    return 0;
}

// linearly interpolate electric field at given points
const int Femocs::interpolate_elfield(int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    check_message(interpolator.get_n_atoms() == 0, "No solution to export!");
    double t0;

    SolutionReader sr(&interpolator);
    start_msg(t0, "=== Interpolating electric field...");
    sr.interpolate(n_points, x, y, z, conf.smoothen_solution*conf.coord_cutoff, 1);
    end_msg(t0);
    if (!skip_calculations) sr.write("output/interpolation_E.movie");

    start_msg(t0, "=== Exporting electric field...");
    sr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    return 0;
}

// linearly interpolate electric potential at given points
const int Femocs::interpolate_phi(int n_points, double* x, double* y, double* z, double* phi, int* flag) {
    check_message(interpolator.get_n_atoms() == 0, "No solution to export!");

    SolutionReader sr(&interpolator);
    sr.interpolate(n_points, x, y, z, conf.smoothen_solution*conf.coord_cutoff, 2, false);
    sr.export_potential(n_points, phi, flag);
    if (!skip_calculations) sr.write("output/interpolation_phi.movie");

    return 0;
}

// parse integer argument of the command from input script
const int Femocs::parse_command(const string& command, int* arg) {
    return conf.read_command(command, arg[0]);
}

// parse double argument of the command from input script
const int Femocs::parse_command(const string& command, double* arg) {
    return conf.read_command(command, arg[0]);
}

} // namespace femocs
