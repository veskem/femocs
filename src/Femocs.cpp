/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"
#include "Coarseners.h"
#include "Macros.h"
#include "Tethex.h"
#include "TetgenMesh.h"

#include <omp.h>
#include <algorithm>

using namespace std;
namespace femocs {

// specify simulation parameters
Femocs::Femocs(const string &conf_file) : skip_calculations(false), fail(false) {
    // Read configuration parameters from configuration file
    conf.read_all(conf_file);

    // Initialise file writing
    MODES.WRITEFILE = conf.n_writefile > 0;

    // Pick the correct verbosity mode flags
    if      (conf.verbose_mode == "mute")    { MODES.QUIET = true;  MODES.VERBOSE = false; }
    else if (conf.verbose_mode == "silent")  { MODES.QUIET = false; MODES.VERBOSE = false; }
    else if (conf.verbose_mode == "verbose") { MODES.QUIET = false; MODES.VERBOSE = true;  }

    // Clear the results from previous run
    if (conf.clear_output) system("rm -rf output");
//    if (MODES.WRITEFILE) system("mkdir -p output");
    system("mkdir -p output");

    start_msg(t0, "======= Femocs started! =======\n");

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_solver  = &ch_solver1;
    prev_ch_solver = NULL;
    end_msg(t0);
}

// delete data and print bye-bye-message
Femocs::~Femocs() {
    start_msg(t0, "======= Femocs finished! =======\n");
}

// Generate boundary nodes for mesh
int Femocs::generate_boundary_nodes(Media& bulk, Media& coarse_surf, Media& vacuum) {
    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract(reader, TYPES.SURFACE);
    end_msg(t0);

    dense_surf.write("output/surface_dense.xyz");

    Coarseners coarseners;
    coarseners.generate(dense_surf, conf.radius, conf.cfactor, conf.latconst);
    coarseners.write("output/coarseners.vtk");

    static bool first_run = true;
    if (first_run) {
        start_msg(t0, "=== Extending surface...");
        if (conf.extended_atoms == "")
            extended_surf = dense_surf.extend(conf.latconst, conf.box_width, coarseners);
        else
            extended_surf = dense_surf.extend(conf.extended_atoms, coarseners);

        end_msg(t0);

        extended_surf.write("output/surface_extended.xyz");
        first_run = false;
    }

    start_msg(t0, "=== Coarsening & smoothing surface...");
    dense_surf.sort_atoms(3, "down");
    coarse_surf = extended_surf;
    coarse_surf += dense_surf;
    coarse_surf = coarse_surf.clean(coarseners);
    coarse_surf.smoothen(conf.radius, conf.surface_smooth_factor, 3.0*conf.coord_cutoff);
    end_msg(t0);

    coarse_surf.write("output/surface_coarse.xyz");

    start_msg(t0, "=== Generating bulk & vacuum...");
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

// Generate bulk and vacuum meshes
int Femocs::generate_meshes(TetgenMesh& bulk_mesh, TetgenMesh& vacuum_mesh) {
    Media bulk, coarse_surf, vacuum;
    fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
    if (fail) return 1;
    
    start_msg(t0, "=== Making big mesh...");
    TetgenMesh big_mesh;
    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality
    fail = big_mesh.generate(bulk, coarse_surf, vacuum, "rnQq" + conf.mesh_quality);
    check_message(fail, "Triangulation failed! Field calculation will be skipped!");
    end_msg(t0);

    start_msg(t0, "=== Making surface faces...");
    big_mesh.generate_appendices();
    end_msg(t0);

    big_mesh.faces.write("output/surface_faces.vtk");

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = big_mesh.mark_mesh(conf.postprocess_marking);
    big_mesh.nodes.write("output/tetmesh_nodes.xyz");
    big_mesh.elems.write("output/tetmesh_elems.vtk");
    check_message(fail, "Mesh marking failed! Field calcualtion will be skipped!");
    end_msg(t0);

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    big_mesh.generate_hexahedra();
    end_msg(t0);

    big_mesh.nodes.write("output/hexmesh_nodes.xyz");
    big_mesh.hexahedra.write("output/hexmesh_elems.vtk");

    start_msg(t0, "=== Separating vacuum & bulk meshes...");
    big_mesh.separate_meshes(bulk_mesh, vacuum_mesh, "rnQ");
    bulk_mesh.group_hexahedra();
    vacuum_mesh.group_hexahedra();
    end_msg(t0);

    start_msg(t0, "=== Cleaning surface faces & atoms...");
    bulk_mesh.faces.clean_sides(reader.sizes);
    dense_surf.clean(bulk_mesh, conf.surface_thichness);
    end_msg(t0);

    bulk_mesh.faces.write("output/surface_faces_clean.vtk");
    dense_surf.write("output/surface_dense_clean.xyz");

    expect(bulk_mesh.nodes.size() > 0, "Zero nodes in bulk mesh!");
    expect(vacuum_mesh.nodes.size() > 0, "Zero nodes in vacuum mesh!");
    expect(bulk_mesh.hexahedra.size() > 0, "Zero elements in bulk mesh!");
    expect(vacuum_mesh.hexahedra.size() > 0, "Zero elements in vacuum mesh!");

    bulk_mesh.hexahedra.write("output/hexmesh_bulk" + conf.message + ".vtk");
    if (MODES.VERBOSE)
        cout << "Bulk:   " << bulk_mesh << "\nVacuum: " << vacuum_mesh << endl;

    return 0;
}

// Solve Laplace equation
int Femocs::solve_laplace(const TetgenMesh& mesh, fch::Laplace<3>& solver) {
    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !solver.import_mesh_directly(mesh.nodes.export_dealii(), mesh.hexahedra.export_dealii());
    check_message(fail, "Importing mesh to Deal.II failed! Field calculation will be skipped!");
    end_msg(t0);
/*
    if (conf.refine_apex) {
        start_msg(t0, "=== Refining mesh in Laplace solver...");
        dealii::Point<3> origin(dense_surf.sizes.xmid, dense_surf.sizes.ymid, dense_surf.sizes.zmax);
        solver.refine_mesh(origin, 7*conf.latconst);
        end_msg(t0);
    }
*/
    start_msg(t0, "=== Initializing Laplace solver...");
    solver.set_applied_efield(conf.neumann);
    solver.setup_system();
    solver.assemble_system();
    end_msg(t0);

    if (MODES.VERBOSE) cout << solver << endl;

    start_msg(t0, "=== Running Laplace solver...");
    solver.solve();
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    vacuum_interpolator.extract_solution(&solver, mesh);
    end_msg(t0);

    vacuum_interpolator.write("output/result_E_phi.xyz");
    vacuum_interpolator.print_statistics();

    return fail;
}

// Solve heat and continuity equations
int Femocs::solve_heat(const TetgenMesh& mesh, fch::Laplace<3>& laplace_solver) {
    if (!conf.heating) return 0;

    start_msg(t0, "=== Initializing rho & T solver...");
    ch_solver->reinitialize(&laplace_solver, prev_ch_solver);
    end_msg(t0);

    start_msg(t0, "=== Importing mesh to rho & T solver...");
    fail = !ch_solver->import_mesh_directly(mesh.nodes.export_dealii(), mesh.hexahedra.export_dealii());
    check_message(fail, "Importing mesh to Deal.II failed! rho & T calculation will be skipped!");
    ch_solver->setup_system();
    end_msg(t0);

    if (MODES.VERBOSE) cout << *(ch_solver) << endl;

    start_msg(t0, "=== Running rho & T solver...\n");
    double t_error = ch_solver->run_specific(conf.t_error, conf.n_newton, false, "", MODES.VERBOSE, 2.0);
    end_msg(t0);

    check_message(t_error > conf.t_error, "Temperature didn't converge, err=" + to_string(t_error)) + "! Using previous solution!";

    start_msg(t0, "=== Extracting T and rho...");
    bulk_interpolator.extract_solution(ch_solver, mesh);
    end_msg(t0);

    bulk_interpolator.write("output/result_rho_T.xyz");

    // Swap current-and-heat-solvers to use solution from current run as a guess in the next one
    static bool odd_run = true;
    if (odd_run) {
        ch_solver = &ch_solver2; prev_ch_solver = &ch_solver1;
    }
    else {
        ch_solver = &ch_solver1; prev_ch_solver = &ch_solver2;
    }
    odd_run = !odd_run;

    return 0;
}

int Femocs::extract_charge(const TetgenMesh& mesh) {
    start_msg(t0, "=== Calculating face charges...");
//    face_charges.calc_interpolated_charges(mesh, conf.E0);  // electric field in the middle of face is interpolated
    face_charges.calc_charges(mesh, conf.E0);             // electric field in the middle of face in directly from solution
    end_msg(t0);

    face_charges.print_statistics(conf.E0 * reader.sizes.xbox * reader.sizes.ybox * face_charges.eps0);
    face_charges.clean(dense_surf.sizes, conf.latconst);
    face_charges.write("output/face_charges.xyz");

    return 0;
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int Femocs::run(const double elfield, const string &message) {
    static unsigned int timestep = 0;  // Counter to measure how many times Femocs has been called
    static bool prev_skip_calculations = true;  // Value of skip_calculations in last call
    double tstart;                     // Variable used to measure the code execution time

    if (!prev_skip_calculations && MODES.WRITEFILE)
        MODES.WRITEFILE = false;

    if ((conf.n_writefile > 0) && (timestep % conf.n_writefile == 0))
        MODES.WRITEFILE = true;

    conf.message = to_string(timestep++);
    write_message("\nFEMOCS: running at timestep " + conf.message);
    conf.message = "_" + string( max(0.0, 5.0 - conf.message.length()), '0' ) + conf.message;

    prev_skip_calculations = skip_calculations;

    check_message(skip_calculations, "Atoms haven't moved significantly, " +
            to_string(reader.rms_distance).substr(0,5) + " < " + to_string(conf.distance_tol).substr(0,5)
            + "! Field calculation will be skipped!");

    skip_calculations = true;
    conf.E0 = elfield; // long-range electric field
    conf.neumann = -10.0 * elfield;  // set minus gradient of solution to equal to E0; also convert V/Angstrom  to  V/nm
    tstart = omp_get_wtime();

    TetgenMesh bulk_mesh;   // FEM mesh in bulk material
    TetgenMesh vacuum_mesh; // FEM mesh in vacuum

    // Generate FEM mesh
    fail = generate_meshes(bulk_mesh, vacuum_mesh);
    if (fail) return 1;

    // Solve Laplace equation on vacuum mesh
    fch::Laplace<3> laplace_solver;
    fail = solve_laplace(vacuum_mesh, laplace_solver);
    if (fail) return 1;

    // Solve heat & continuity equation on bulk mesh
    fail = solve_heat(bulk_mesh, laplace_solver);
    if (fail) return 1;

    // Extract face charges
    fail = extract_charge(bulk_mesh);
    if (fail) return 1;

    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.distance_tol);
    end_msg(t0);

    write_message("\nFEMOCS: total execution time " + to_string(omp_get_wtime() - tstart));
    skip_calculations = false;

    return 0;
}

// import atoms from file
int Femocs::import_atoms(const string& file_name) {
    clear_log();
    string file_type, fname;

    if (file_name == "") fname = conf.atom_file;
    else fname = file_name;

    file_type = get_file_type(fname);
    expect(file_type == "ckx" || file_type == "xyz", "Unknown file type: " + file_type);

    start_msg(t0, "=== Importing atoms...");
    reader.import_file(fname);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.size() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        if (file_type == "xyz") {
            start_msg(t0, "=== Generating list of close neighbours...");
            reader.calc_nborlist(conf.nnn, conf.coord_cutoff);
            end_msg(t0);

            start_msg(t0, "=== Coordination & cluster analysis...");
            reader.calc_coordinations();
            if (conf.cluster_anal) reader.calc_clusters();
            end_msg(t0);
            reader.check_clusters();

            start_msg(t0, "=== Extracting atom types...");
            reader.extract_types(conf.nnn, conf.latconst);
            end_msg(t0);

        } else {
            start_msg(t0, "=== Calculating coords from atom types...");
            reader.calc_coordinations(conf.nnn);
            end_msg(t0);
        }
        reader.write("output/atomreader.xyz");
    }

    return 0;
}

// import atoms from PARCAS
int Femocs::import_atoms(const int n_atoms, const double* coordinates, const double* box, const int* nborlist) {
    clear_log();

    start_msg(t0, "=== Importing atoms...");
    reader.import_parcas(n_atoms, coordinates, box);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.size() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Generating list of close neighbours...");
        reader.calc_nborlist(conf.nnn, conf.coord_cutoff, nborlist);
        end_msg(t0);

        start_msg(t0, "=== Coordination & cluster analysis...");
        reader.calc_coordinations();
        if (conf.cluster_anal) reader.calc_clusters();
        end_msg(t0);
        reader.check_clusters();

        start_msg(t0, "=== Extracting atom types...");
        reader.extract_types(conf.nnn, conf.latconst);
        end_msg(t0);

        reader.write("output/atomreader.xyz");
    }

    return 0;
}

// import coordinates and types of atoms
int Femocs::import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    clear_log();

    start_msg(t0, "=== Importing atoms...");
    reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);
    if (MODES.VERBOSE) cout << "#input atoms: " << reader.size() << endl;

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating coords from atom types...");
        reader.calc_coordinations(conf.nnn);
        end_msg(t0);

        if (conf.cluster_anal) reader.calc_clusters();
        reader.write("output/atomreader.xyz");
    }

    return 0;
}

// calculate and export electric field on imported atom coordinates
int Femocs::export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    check_message(vacuum_interpolator.size() == 0, "No field to export!");

    if (!skip_calculations) {
        start_msg(t0, "=== Interpolating E and phi...");
        fields.interpolate(dense_surf, conf.use_histclean * conf.coord_cutoff);
        end_msg(t0);

        fields.write("output/elfields.movie");
        fields.print_statistics();
    }

    start_msg(t0, "=== Exporting electric field...");
    fields.export_solution(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);

    return 0;
}

// calculate and export temperatures on imported atom coordinates
int Femocs::export_temperature(const int n_atoms, double* T) {
    if (!conf.heating) return 0;
    check_message(bulk_interpolator.size() == 0, "No temperature to export!");

    if (!skip_calculations) {
        start_msg(t0, "=== Interpolating T and rho...");
        temperatures.interpolate(reader, 0);
        end_msg(t0);

        temperatures.write("output/interpolation_bulk.movie");
        temperatures.print_statistics();
    }

    start_msg(t0, "=== Exporting temperatures...");
    temperatures.export_temperature(n_atoms, T);
    end_msg(t0);

    return 0;
}

// calculate and export charges & forces on imported atom coordinates
int Femocs::export_charge_and_force(const int n_atoms, double* xq) {
    check_message(fields.size() == 0 || face_charges.size() == 0, "No force to export!");

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating atomic forces...");
        forces.calc_forces(fields, face_charges, conf.use_histclean*conf.coord_cutoff, conf.charge_smooth_factor);
        end_msg(t0);

        forces.write("output/forces.movie");
        forces.print_statistics(conf.E0 * reader.sizes.xbox * reader.sizes.ybox * face_charges.eps0);
    }

    start_msg(t0, "=== Exporting atomic forces...");
    forces.export_force(n_atoms, xq);
    end_msg(t0);

    return 0;
}

// linearly interpolate electric field at given points
int Femocs::interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    check_message(vacuum_interpolator.size() == 0, "No solution to export!");

    FieldReader fr(&vacuum_interpolator);
    start_msg(t0, "=== Interpolating electric field...");
    fr.interpolate(n_points, x, y, z, conf.use_histclean * conf.coord_cutoff, 1);
    end_msg(t0);
    if (!skip_calculations) fr.write("output/interpolation_E.movie");

    start_msg(t0, "=== Exporting electric field...");
    fr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    return 0;
}

// linearly interpolate electric potential at given points
int Femocs::interpolate_phi(const int n_points, const double* x, const double* y, const double* z, double* phi, int* flag) {
    check_message(vacuum_interpolator.size() == 0, "No solution to export!");

    FieldReader fr(&vacuum_interpolator);
    fr.interpolate(n_points, x, y, z, conf.use_histclean * conf.coord_cutoff, 2, false);
    fr.export_potential(n_points, phi, flag);
    if (!skip_calculations) fr.write("output/interpolation_phi.movie");

    return 0;
}

// parse integer argument of the command from input script
int Femocs::parse_command(const string& command, int* arg) {
    return conf.read_command(command, arg[0]);
}

// parse double argument of the command from input script
int Femocs::parse_command(const string& command, double* arg) {
    return conf.read_command(command, arg[0]);
}

} // namespace femocs
