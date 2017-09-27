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
#include "VoronoiMesh.h"

#include <omp.h>
#include <algorithm>
#include <sstream>

using namespace std;
namespace femocs {

// specify simulation parameters
Femocs::Femocs(const string &conf_file) : skip_calculations(false), fail(false), timestep(-1) {
    static bool first_call = true;

    // Read configuration parameters from configuration file
    conf.read_all(conf_file);

    // Initialise file writing
    MODES.WRITEFILE = conf.n_writefile > 0;

    // Pick the correct verbosity mode flags
    if      (conf.verbose_mode == "mute")    { MODES.MUTE = true;  MODES.VERBOSE = false; }
    else if (conf.verbose_mode == "silent")  { MODES.MUTE = false; MODES.VERBOSE = false; }
    else if (conf.verbose_mode == "verbose") { MODES.MUTE = false; MODES.VERBOSE = true;  }

    // Clear the results from previous run
    if (first_call && conf.clear_output) fail = system("rm -rf out");
    fail = system("mkdir -p out");

    start_msg(t0, "======= Femocs started! =======\n");

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_solver  = &ch_solver1;

    ch_transient_solver.set_physical_quantities(&phys_quantities);

    prev_ch_solver = NULL;
    end_msg(t0);

    first_call = false;
}

// delete data and print bye-bye-message
Femocs::~Femocs() {
    start_msg(t0, "======= Femocs finished! =======\n");
}

// Write all the available data to file for debugging purposes
int Femocs::force_output() {
    if (conf.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    fem_mesh.hexahedra.write("out/hexmesh_err.vtk");
    fem_mesh.elems.write("out/tetmesh_err.vtk");
    fem_mesh.faces.write("out/trimesh_err.vtk");

    vacuum_interpolator.write("out/result_E_phi_vacuum_err.xyz");
    vacuum_interpolator.write("out/result_E_phi_vacuum_err.vtk");
    surface_interpolator.write("out/result_E_phi_surface_err.xyz");
    surface_interpolator.write("out/result_E_phi_surface_err.vtk");

    if (bulk_interpolator.size() > 0) {
        if (conf.heating_mode == "transient" || conf.heating_mode == "stationary")
            bulk_interpolator.write("out/result_J_T_err.xyz");
        if (conf.heating_mode == "transient") {
            ch_transient_solver.output_results_current("out/result_J_err.vtk");
            ch_transient_solver.output_results_heating("out/result_T_err.vtk");
        }
        if (conf.heating_mode == "stationary") {
            ch_solver->output_results("out/result_J_T_err.vtk");
        }
    }

    return 0;
}

// Interpolate the solution on the x-z plane in the middle of simulation box
void Femocs::write_slice(const string& file_name) {
    int writefile_save = MODES.WRITEFILE;
    MODES.WRITEFILE = true;

    const int nx = 300;  // number of points in x-direction
    const int nz = 300;  // number of points in z-direction
    const double eps = 1e-5 * conf.latconst;

    const double xmax = reader.sizes.xmid;
    const double xmin = xmax - 3*conf.radius;
    const double zmin = reader.sizes.zmin;
    const double zmax = zmin + 3*conf.radius;
    const double dx = (xmax - xmin) / (nx-1);
    const double dz = (zmax - zmin) / (nz-1);

    Medium medium(nx * nz);

    for (double x = xmin; x < xmax + eps; x += dx)
        for (double z = zmin; z < zmax + eps; z += dz)
            medium.append( Point3(x, reader.sizes.ymid, z) );

    FieldReader fr(&vacuum_interpolator);
    fr.interpolate(medium, conf.use_histclean * conf.coordination_cutoff);
    fr.write(file_name);

    MODES.WRITEFILE = writefile_save;
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int Femocs::run(const double elfield, const string &message) {
    stringstream stream; stream << fixed << setprecision(3);

    // convert message to integer time step
    int tstep;
    stream << message;
    if (!stream >> tstep)
        tstep = -1;

    stream.str("");
    stream << "Atoms haven't moved significantly, " << reader.rms_distance
        << " < " << conf.distance_tol << "! Field calculation will be skipped!";
    check_return(reinit(tstep), stream.str());

    double tstart = omp_get_wtime();
    skip_calculations = true;

    check_return(generate_meshes(), "Mesh generation failed!");

    // Solve Laplace equation on vacuum mesh
    if (solve_laplace(elfield)) {
        force_output();
        check_return(true, "Solving Laplace equation failed!");
    }

    // Solve heat & continuity equation on bulk mesh
    if (solve_heat(conf.t_ambient)) {
        force_output();
        check_return(true, "Solving heat & continuity equation failed!");
    }

    finalize();

    stream.str("");
    stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());

    return 0;
}

// Determine whether atoms have moved significantly and whether to enable file writing
int Femocs::reinit(const int tstep) {
    static bool prev_skip_calculations = true;  // Value of skip_calculations in last call
    if (tstep >= 0)
        timestep = tstep;
    else
        ++timestep;

    if (!prev_skip_calculations && MODES.WRITEFILE)
        MODES.WRITEFILE = false;

    if ((conf.n_writefile > 0) && (timestep % conf.n_writefile == 0))
        MODES.WRITEFILE = true;

    conf.message = to_string(timestep);
    write_silent_msg("Running at timestep " + conf.message);
    conf.message = "_" + string( max(0.0, 6.0 - conf.message.length()), '0' ) + conf.message;

    prev_skip_calculations = skip_calculations;
    return skip_calculations;
}

// Store the imported atom coordinates and set flag that enables exporters
int Femocs::finalize() {
    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.distance_tol);
    end_msg(t0);
    skip_calculations = false;
    return 0;
}

// Generate boundary nodes for mesh
int Femocs::generate_boundary_nodes(Media& bulk, Media& coarse_surf, Media& vacuum) {
    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract(reader, TYPES.SURFACE);
    dense_surf.sort_atoms(3, "down");
    end_msg(t0);

    dense_surf.write("out/surface_dense.xyz");

    if (conf.surface_cleaner == "voronois") {
        start_msg(t0, "=== Cleaning surface with Voronoi cells...");
//        const int err_code = dense_surf.voronoi_clean(areas, conf.radius, conf.latconst, conf.mesh_quality + "a10");
        const int err_code = dense_surf.voronoi_clean(areas, conf.radius, conf.latconst, "1.4");
        check_return(err_code, "Making voronoi cells failed with error code " + to_string(err_code));
        end_msg(t0);
        dense_surf.write("out/surface_dense_clean.xyz");
    }

    coarseners.generate(dense_surf, conf.radius, conf.cfactor, conf.latconst);
    coarseners.write("out/coarseners.vtk");

    static bool first_run = true;
    if (first_run) {
        start_msg(t0, "=== Extending surface...");
        if (conf.extended_atoms == "")
            extended_surf = dense_surf.extend(conf.latconst, conf.box_width, coarseners);
        else
            extended_surf = dense_surf.extend(conf.extended_atoms, coarseners);

        end_msg(t0);

        extended_surf.write("out/surface_extended.xyz");
        first_run = false;
    }

    start_msg(t0, "=== Coarsening & smoothing surface...");
    coarse_surf = extended_surf;
    coarse_surf += dense_surf;
    coarse_surf = coarse_surf.clean(coarseners);

    coarse_surf.smoothen(conf.radius, conf.surface_smooth_factor, 3.0*conf.coordination_cutoff);
//    coarse_surf.smoothen(conf, 3.0 * conf.surface_smooth_factor * conf.coordination_cutoff);
    end_msg(t0);

    coarse_surf.write("out/surface_coarse.xyz");

    start_msg(t0, "=== Generating bulk & vacuum...");
    coarse_surf.calc_statistics();  // calculate zmin and zmax for surface
    vacuum.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin + conf.box_height * coarse_surf.sizes.zbox);
    bulk.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin - conf.bulk_height * conf.latconst);
    reader.resize_box(coarse_surf.sizes.xmin, coarse_surf.sizes.xmax, 
        coarse_surf.sizes.ymin, coarse_surf.sizes.ymax,
        bulk.sizes.zmin, vacuum.sizes.zmax);
    end_msg(t0);
    
    bulk.write("out/bulk.xyz");
    vacuum.write("out/vacuum.xyz");

    return 0;
}

// Generate bulk and vacuum meshes
int Femocs::generate_meshes() {
    fem_mesh.clear();

    Media bulk, coarse_surf, vacuum;
    fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
    if (fail) return 1;

    start_msg(t0, "=== Making big mesh...");
    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality, a - element volume,
    // F - suppress output of faces and edges, B - suppress output of boundary info
    string command = "rnQFBq" + conf.mesh_quality;
    if (conf.element_volume != "") command += "a" + conf.element_volume;
    int err_code = fem_mesh.generate(bulk, coarse_surf, vacuum, command);
    check_return(err_code, "Triangulation failed with error code " + to_string(err_code));
    end_msg(t0);

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = fem_mesh.mark_mesh();
    fem_mesh.nodes.write("out/tetmesh_nodes.xyz");
    fem_mesh.elems.write("out/tetmesh.vtk");
    check_return(fail, "Mesh marking failed!");
    end_msg(t0);

    start_msg(t0, "=== Generating surface faces...");
    fail = fem_mesh.generate_surface(reader.sizes, "rnQB");
    end_msg(t0);
    fem_mesh.faces.write("out/trimesh.vtk");
    check_return(fail, "Generation of surface faces failed!");

    if (conf.smooth_algorithm != "none" && conf.smooth_steps > 0) {
        start_msg(t0, "=== Smoothing triangles...");
        fem_mesh.smoothen_tris(conf.smooth_steps, conf.smooth_lambda, conf.smooth_mu, conf.smooth_algorithm);
        end_msg(t0);
        fem_mesh.faces.write("out/trimesh_smooth.vtk");
    }

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    fem_mesh.generate_hexahedra();
    end_msg(t0);

    fem_mesh.nodes.write("out/hexmesh_nodes.xyz");
    fem_mesh.quads.write("out/quadmesh.vtk");
    fem_mesh.hexahedra.write("out/hexmesh.vtk");
    fem_mesh.write_separate("out/hexmesh_bulk" + conf.message + ".vtk", TYPES.BULK);
    stringstream ss; ss << fem_mesh;
    write_verbose_msg(ss.str());

    if (conf.surface_cleaner == "faces") {
        start_msg(t0, "=== Cleaning surface with triangles...");
        dense_surf.faces_clean(fem_mesh, conf.surface_thichness);
        end_msg(t0);

        dense_surf.write("out/surface_dense_clean.xyz");
    }

    return 0;
}

// Solve Laplace equation
int Femocs::solve_laplace(const double E0) {
    conf.E0 = E0;       // long-range electric field
    conf.neumann = -E0; // set minus gradient of solution to equal to E0

    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(E0, conf.field_tolerance_min, conf.field_tolerance_max, conf.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !laplace_solver.import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_vacuum());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Initializing Laplace solver...");
    laplace_solver.set_applied_efield(conf.neumann);
    laplace_solver.setup_system();
    laplace_solver.assemble_system();
    end_msg(t0);

    stringstream ss; ss << laplace_solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Laplace solver...");
    laplace_solver.solve(conf.n_phi, conf.phi_error, true, conf.ssor_param);
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    fail = vacuum_interpolator.extract_solution(&laplace_solver);
    fail |= surface_interpolator.extract_solution(&laplace_solver);
    end_msg(t0);

    check_return(fields.check_limits(vacuum_interpolator.get_solutions()), "Vacuum field is over-enhanced!");
    check_return(fields.check_limits(surface_interpolator.get_solutions()), "Surface field is over-enhanced!");

    vacuum_interpolator.write("out/result_E_phi_vacuum.xyz");
    vacuum_interpolator.write("out/result_E_phi_vacuum.vtk");
    surface_interpolator.write("out/result_E_phi_surface.xyz");
    surface_interpolator.write("out/result_E_phi_surface.vtk");

    return fail;
}

// Pick a method to solve heat & continuity equations
int Femocs::solve_heat(const double T_ambient) {
    if (conf.heating_mode == "stationary")
        return solve_stationary_heat(T_ambient);

    else if (conf.heating_mode == "transient")
        for (int i = 0; i < conf.transient_steps; ++i) {
            if (solve_transient_heat(T_ambient))
                return 1;
        }

    return 0;
}

// Solve steady-state heat and continuity equations
int Femocs::solve_stationary_heat(const double T_ambient) {
    start_msg(t0, "=== Initializing stationary J & T solver...");
    ch_solver->reinitialize(&laplace_solver, prev_ch_solver);
    end_msg(t0);

    stringstream ss; ss << *(ch_solver);
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Importing mesh to J & T solver...");
    fail = !ch_solver->import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_bulk());
    check_return(fail, "Importing mesh to Deal.II failed!");
    ch_solver->setup_system();
    end_msg(t0);

    start_msg(t0, "=== Transfering elfield to J & T solver...");
    FieldReader fr(&vacuum_interpolator);
    fr.transfer_elfield(ch_solver, conf.coordination_cutoff, conf.use_histclean);
    end_msg(t0);

    fr.write("out/surface_field.xyz");

    start_msg(t0, "=== Running J & T solver...\n");
    ch_solver->set_ambient_temperature(T_ambient);
    double t_error = ch_solver->run_specific(conf.t_error, conf.n_newton, false, "", MODES.VERBOSE, 2, 400, true);
    end_msg(t0);

    ch_solver->output_results("out/result_J_T.vtk");

    check_return(t_error > conf.t_error, "Temperature didn't converge, err=" + to_string(t_error));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator.write("out/result_J_T.xyz");

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

// Solve transient heat and continuity equations
int Femocs::solve_transient_heat(const double T_ambient) {
    static bool first_call = true;
    const double time_unit = 1e-12; // == picosec
    double delta_time = 0.01 * time_unit;
    const int n_time_steps = conf.transient_time / delta_time;

    require(n_time_steps > 0, "Too small transient_time: " + to_string(conf.transient_time));
    delta_time = conf.transient_time / n_time_steps;

    start_msg(t0, "=== Importing mesh to transient J & T solver...");
    fail = !ch_transient_solver.import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_bulk());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Transfering elfield to J & T solver...");
    FieldReader field_reader(&vacuum_interpolator);
    field_reader.transfer_elfield(ch_transient_solver, conf.coordination_cutoff, conf.use_histclean);
    end_msg(t0);
    field_reader.write("out/surface_field.xyz");

    start_msg(t0, "=== Interpolating J & T on face centroids...");
    HeatReader heat_reader(&bulk_interpolator);
    heat_reader.interpolate(ch_transient_solver, T_ambient);
    end_msg(t0);
    heat_reader.write("out/surface_temperature.xyz");

    start_msg(t0, "=== Calculating field emission...");
    EmissionReader emission(&vacuum_interpolator);
    emission.transfer_emission(ch_transient_solver, field_reader, conf.work_function, heat_reader);
    end_msg(t0);
    emission.write("out/surface_emission.xyz");

    if (first_call) {
        start_msg(t0, "=== Setup transient J & T solver...");
        ch_transient_solver.setup_current_system();
        ch_transient_solver.setup_heating_system();
        end_msg(t0);
    }

    start_msg(t0, "=== Calculating current density...");
    ch_transient_solver.assemble_current_system();           // assemble matrix for current density equation; current == electric current
    unsigned int ccg = ch_transient_solver.solve_current();  // ccg == number of current calculation (CG) iterations
    end_msg(t0);

    start_msg(t0, "=== Calculating temperature distribution...\n");
    ch_transient_solver.set_timestep(delta_time);

    for (int i = 0; i < n_time_steps; ++i) {
        ch_transient_solver.assemble_heating_system_euler_implicit();
        //ch_transient_solver.assemble_heating_system_crank_nicolson();

        unsigned int hcg = ch_transient_solver.solve_heat(); // hcg == number of temperature calculation (CG) iterations
        if (MODES.VERBOSE) {
            double max_T = ch_transient_solver.get_max_temperature();
            printf("  t=%5.3fps; ccg=%2d; hcg=%2d; Tmax=%6.2f\n", i * delta_time / time_unit, ccg, hcg, max_T);
        }
    }
    end_msg(t0);

    ch_transient_solver.output_results_current("out/result_J.vtk");
    ch_transient_solver.output_results_heating("out/result_T.vtk");

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.extract_solution(&ch_transient_solver);
    end_msg(t0);
    bulk_interpolator.write("out/result_J_T.movie");

    first_call = false;
    return 0;
}

// Generate artificial nanotip
int Femocs::generate_nanotip(const double height, const double radius, const double resolution) {
    clear_log();

    double res = conf.latconst;
    if (resolution > 0) res = resolution;

    double r = conf.radius - res;
    if (radius > 0) r = radius;

    start_msg(t0, "=== Generating nanotip...");
    reader.generate_nanotip(height, r, res);
    reader.calc_coordinations(conf.nnn);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    reader.write("out/atomreader.ckx");
    return 0;
}

// import atoms from a file
int Femocs::import_atoms(const string& file_name, const int add_noise) {
    clear_log();
    string file_type, fname;

    if (file_name == "") fname = conf.atom_file;
    else fname = file_name;

    file_type = get_file_type(fname);
    require(file_type == "ckx" || file_type == "xyz", "Unknown file type: " + file_type);

    start_msg(t0, "=== Importing atoms...");
    reader.import_file(fname, add_noise);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        if (file_type == "xyz") {
            start_msg(t0, "=== Performing coordination analysis...");
            reader.calc_coordinations(conf.nnn, conf.coordination_cutoff);
            end_msg(t0);

            if (conf.cluster_anal) {
                start_msg(t0, "=== Performing cluster analysis...");
                if (conf.cluster_cutoff <= 0) reader.calc_clusters();
                else reader.calc_clusters(conf.nnn, conf.cluster_cutoff);
                end_msg(t0);
                reader.check_clusters(1);
            }

            start_msg(t0, "=== Extracting atom types...");
            reader.extract_types(conf.nnn, conf.latconst);
            end_msg(t0);

        } else {
            start_msg(t0, "=== Calculating coords from atom types...");
            reader.calc_coordinations(conf.nnn);
            end_msg(t0);
        }
    }

    reader.write("out/atomreader.ckx");
    return 0;
}

// import atoms from PARCAS
int Femocs::import_atoms(const int n_atoms, const double* coordinates, const double* box, const int* nborlist) {
    clear_log();

    start_msg(t0, "=== Importing atoms...");
    reader.import_parcas(n_atoms, coordinates, box);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Performing coordination analysis...");
        reader.calc_coordinations(conf.nnn, conf.coordination_cutoff, nborlist);
        end_msg(t0);

        if (conf.cluster_anal) {
            start_msg(t0, "=== Performing cluster analysis...");
            if (conf.cluster_cutoff <= 0) reader.calc_clusters();
            else reader.calc_clusters(conf.nnn, conf.cluster_cutoff, nborlist);
            end_msg(t0);
            reader.check_clusters(1);
        }

        start_msg(t0, "=== Extracting atom types...");
        reader.extract_types(conf.nnn, conf.latconst);
        end_msg(t0);
    }

    reader.write("out/atomreader.ckx");
    return 0;
}

// import coordinates and types of atoms
int Femocs::import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    clear_log();
    conf.surface_cleaner = "none"; // disable the surface cleaner for atoms with known types

    start_msg(t0, "=== Importing atoms...");
    reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.distance_tol) < conf.distance_tol;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating coordinations from atom types...");
        reader.calc_coordinations(conf.nnn);
        end_msg(t0);
    }

    reader.write("out/atomreader.ckx");
    return 0;
}

// export the atom types as seen by FEMOCS
int Femocs::export_atom_types(const int n_atoms, int* types) {
    const int export_size = min(n_atoms, reader.size());
    for (int i = 0; i < export_size; ++i)
        types[i] = reader.get_marker(i);

    return reader.check_clusters(0);
}

// calculate and export electric field on imported atom coordinates
int Femocs::export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm) {
    if (n_atoms < 0) return 0;
    check_return(vacuum_interpolator.size() == 0, "No field to export!");

    fail = false;

    if (skip_calculations)
        write_silent_msg("Using previous electric field!");
    else {
        start_msg(t0, "=== Interpolating E and phi...");
        fields.interpolate2D(dense_surf, 0, true);
        fail = fields.clean(conf.coordination_cutoff, conf.use_histclean);
        end_msg(t0);

        fields.write("out/fields.movie");
        fail |= fields.check_limits();
    }

    start_msg(t0, "=== Exporting electric field...");
    fields.export_solution(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);

    return fail;
}

// calculate and export temperatures on imported atom coordinates
int Femocs::export_temperature(const int n_atoms, double* T) {
    if (n_atoms < 0 || conf.heating_mode == "none") return 0;
    check_return(bulk_interpolator.size() == 0, "No temperature to export!");

    if (skip_calculations)
        write_silent_msg("Using previous temperature!");
    else {
        start_msg(t0, "=== Interpolating J & T...");
        temperatures.interpolate(reader);
        end_msg(t0);

        temperatures.write("out/interpolation_bulk.movie");
        temperatures.print_statistics();
    }

    start_msg(t0, "=== Exporting temperature...");
    temperatures.export_temperature(n_atoms, T);
    end_msg(t0);

    return 0;
}

// calculate and export charges & forces on imported atom coordinates
int Femocs::export_charge_and_force(const int n_atoms, double* xq) {
    if (n_atoms < 0) return 0;
    check_return(fields.size() == 0, "No charge & force to export!");

    if (skip_calculations)
        write_silent_msg("Using previous charge & force!");
    else {
        // analytical total charge without epsilon0 (will be added in ChargeReader)
        const double tot_charge = conf.E0 * reader.sizes.xbox * reader.sizes.ybox;

        ChargeReader face_charges(&vacuum_interpolator); // charges on surface triangles
        face_charges.set_check_params(tot_charge, conf.charge_tolerance_min, conf.charge_tolerance_max);

        start_msg(t0, "=== Calculating face charges...");
        face_charges.calc_interpolated_charges(fem_mesh, conf.E0);
//        face_charges.calc_charges(fem_mesh, conf.E0);
        end_msg(t0);
//        check_return(face_charges.check_limits(), "Face charges are not conserved!");

        face_charges.clean(dense_surf.sizes, conf.latconst);
        face_charges.write("out/face_charges.xyz");

        start_msg(t0, "=== Calculating charges & forces...");
        forces.calc_forces(fields, face_charges, conf.use_histclean*conf.coordination_cutoff,
                conf.charge_smooth_factor);
//        forces.calc_forces_vol2(vacuum_mesh, fields, face_charges, surface_interpolator,
//                conf.use_histclean*conf.coordination_cutoff, conf.charge_smooth_factor);
        end_msg(t0);
        face_charges.check_limits(forces.get_interpolations());
        forces.write("out/forces_before.xyz");

        start_msg(t0, "=== Calculating Voronoi charges & forces...");
//        forces.recalc_forces(fields, areas);
//        forces.calc_voronoi_charges(conf.radius, conf.latconst, "1.1");
        forces.calc_surface_voronoi_charges(fem_mesh.elems, fields, conf.radius, conf.latconst, "1.2");
//        forces.calc_kmc_voronoi_charges(reader, fem_mesh.elems, fields, conf.radius, conf.latconst, "10.0");
        end_msg(t0);

        forces.write("out/forces.movie");
        face_charges.check_limits(forces.get_interpolations());
    }

    start_msg(t0, "=== Exporting atomic charges & forces...");
    forces.export_force(n_atoms, xq);
    end_msg(t0);

    return 0;
}

// linearly interpolate electric field at given points
int Femocs::interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    if (n_points <= 0) return 0;
    check_return(surface_interpolator.size() == 0, "No solution to export!");

    FieldReader fr(&surface_interpolator);
    start_msg(t0, "=== Interpolating & exporting surface elfield...");
    fr.interpolate2D(n_points, x, y, z, 1, false);
    fail = fr.clean(conf.coordination_cutoff, conf.use_histclean);
    fr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    fr.write("out/interpolation_E_surf.movie");
    return fields.check_limits(fr.get_interpolations()) || fail;
}

// linearly interpolate electric field at given points
int Femocs::interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator.size() == 0, "No electric field to export!");

    FieldReader fr(&vacuum_interpolator);
    start_msg(t0, "=== Interpolating & exporting elfield...");
    fr.interpolate(n_points, x, y, z, 1, false);
    fail = fr.clean(conf.coordination_cutoff, conf.use_histclean);
    fr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    fr.write("out/interpolation_E.movie");
    return fields.check_limits(fr.get_interpolations()) || fail;
}

// linearly interpolate electric potential at given points
int Femocs::interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
        double* phi, int* flag) {

    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator.size() == 0, "No electric potential to export!");

    FieldReader fr(&vacuum_interpolator);
    fr.interpolate(n_points, x, y, z, 2, false);
    fail = fr.clean(conf.coordination_cutoff, conf.use_histclean);
    fr.export_potential(n_points, phi, flag);

    return fail;
}

// parse integer argument of the command from input script
int Femocs::parse_command(const string& command, int* arg) {
    return conf.read_command(command, arg[0]);
}

// parse double argument of the command from input script
int Femocs::parse_command(const string& command, double* arg) {
    return conf.read_command(command, arg[0]);
}

// parse string argument of the command from input script
int Femocs::parse_command(const string& command, string& arg) {
    return conf.read_command(command, arg);
}

// parse char array argument of the command from input script
int Femocs::parse_command(const string& command, char* arg) {
    string string_arg;
    bool fail = conf.read_command(command, string_arg);
    if (!fail) string_arg.copy(arg, string_arg.length());
    return fail;
}

} // namespace femocs
