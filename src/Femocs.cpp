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
Femocs::Femocs(const string &conf_file) : skip_calculations(false), fail(false),
        timestep(-1), last_full_timestep(0) {
    static bool first_call = true;

    // Read configuration parameters from configuration file
    conf.read_all(conf_file);

    // Initialise file writing
    MODES.WRITEFILE = conf.behaviour.n_writefile > 0;

    // Pick the correct verbosity mode flags
    if      (conf.behaviour.verbosity == "mute")    { MODES.MUTE = true;  MODES.VERBOSE = false; }
    else if (conf.behaviour.verbosity == "silent")  { MODES.MUTE = false; MODES.VERBOSE = false; }
    else if (conf.behaviour.verbosity == "verbose") { MODES.MUTE = false; MODES.VERBOSE = true;  }

    // Clear the results from previous run
    if (first_call && conf.run.output_cleaner) fail = system("rm -rf out");
    fail = system("mkdir -p out");
    first_call = false;

    start_msg(t0, "======= Femocs started! =======\n");

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    phys_quantities.initialize_with_hc_data();
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_solver = &ch_solver1;
    prev_ch_solver = NULL;
    ch_transient_solver.set_physical_quantities(&phys_quantities);
    end_msg(t0);

    // Pick the ranks and types of interpolators
    if (conf.behaviour.interpolation_rank == 2) {
        vacuum_surface_interpolator = new QuadTriInterpolator(&fem_mesh);
        bulk_surface_interpolator = new QuadTriInterpolator(&fem_mesh);
        vacuum_interpolator = new QuadTetInterpolator(&fem_mesh);
        bulk_interpolator = new QuadTetInterpolator(&fem_mesh);
    } else {
        vacuum_surface_interpolator = new LinTriInterpolator(&fem_mesh);
        bulk_surface_interpolator = new LinTriInterpolator(&fem_mesh);
        vacuum_interpolator = new LinTetInterpolator(&fem_mesh);
        bulk_interpolator = new LinTetInterpolator(&fem_mesh);
    }

    temperatures.set_interpolator(bulk_interpolator);
    fields.set_interpolator(vacuum_interpolator, vacuum_surface_interpolator);
    forces.set_interpolator(vacuum_interpolator, vacuum_surface_interpolator);
}

// delete data and print bye-bye-message
Femocs::~Femocs() {
    delete vacuum_surface_interpolator;
    delete bulk_surface_interpolator;
    delete vacuum_interpolator;
    delete bulk_interpolator;

    start_msg(t0, "======= Femocs finished! =======\n");
}

// Write all the available data to file for debugging purposes
int Femocs::force_output() {
    if (conf.behaviour.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    fem_mesh.hexahedra.write("out/hexmesh_err.vtk");
    fem_mesh.elems.write("out/tetmesh_err.vtk");
    fem_mesh.faces.write("out/trimesh_err.vtk");

    vacuum_interpolator->write("out/result_E_phi_vacuum_err.xyz");
    vacuum_interpolator->write("out/result_E_phi_vacuum_err.vtk");
    vacuum_surface_interpolator->write("out/result_E_phi_surface_err.xyz");
    vacuum_surface_interpolator->write("out/result_E_phi_surface_err.vtk");

    if (bulk_interpolator->size() > 0) {
        if (conf.heating.mode == "transient" || conf.heating.mode == "stationary")
            bulk_interpolator->write("out/result_J_T_err.xyz");
        if (conf.heating.mode == "transient") {
            ch_transient_solver.output_results_current("out/result_J_err.vtk");
            ch_transient_solver.output_results_heating("out/result_T_err.vtk");
        }
        if (conf.heating.mode == "stationary") {
            ch_solver->output_results("out/result_J_T_err.vtk");
        }
    }

    return 0;
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int Femocs::run(const double elfield, const string &timestep) {
    stringstream parser, stream;
    stream << fixed << setprecision(3);

    // convert message to integer time step
    int tstep;
    parser << timestep;
    parser >> tstep;

    stream.str("");
    stream << "Atoms haven't moved significantly, " << reader.rms_distance
        << " < " << conf.tolerance.distance << "! Field calculation will be skipped!";
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
    if (solve_heat(conf.heating.t_ambient)) {
        force_output();
        check_return(true, "Solving heat & continuity equation failed!");
    }

    finalize();

    stream.str("");
    stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());

    return 0;
}

int Femocs::run() {
    return run(conf.laplace.E0, "");
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

    if ((conf.behaviour.n_writefile > 0) && (timestep % conf.behaviour.n_writefile == 0))
        MODES.WRITEFILE = true;

    atom2face.clear();

    timestep_string = to_string(timestep);
    write_silent_msg("Running at timestep " + timestep_string);
    timestep_string = "_" + string( max(0.0, 6.0 - timestep_string.length()), '0' ) + timestep_string;

    prev_skip_calculations = skip_calculations;
    return skip_calculations;
}

// Store the imported atom coordinates and set flag that enables exporters
int Femocs::finalize() {
    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.tolerance.distance);
    end_msg(t0);
    skip_calculations = false;
    last_full_timestep = timestep;
    return 0;
}

// Generate boundary nodes for mesh
int Femocs::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum) {
    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract(reader, TYPES.SURFACE);
    end_msg(t0);

    dense_surf.write("out/surface_dense.xyz");

    coarseners.generate(dense_surf, conf.geometry.radius, conf.cfactor, conf.geometry.latconst);
    coarseners.write("out/coarseners.vtk");

    static bool first_run = true;
    if (first_run) {
        start_msg(t0, "=== Extending surface...");
        if (conf.path.extended_atoms == "")
            dense_surf.extend(extended_surf, coarseners, conf.geometry.latconst, conf.geometry.box_width);
        else
            extended_surf = dense_surf.extend(conf.path.extended_atoms, coarseners);
        end_msg(t0);
        first_run = false;
    }

    start_msg(t0, "=== Coarsening & smoothing surface...");
    coarse_surf = extended_surf;
//    coarse_surf += dense_surf;
    coarse_surf += dense_surf.clean_roi(coarseners);
    coarse_surf = coarse_surf.clean(coarseners);
    coarse_surf.smoothen(conf.geometry.radius, conf.smoothing.beta_atoms, 3.0*conf.geometry.coordination_cutoff);
    end_msg(t0);

    coarse_surf.write("out/surface_coarse.xyz");

    start_msg(t0, "=== Generating bulk & vacuum corners...");
    coarse_surf.calc_statistics();  // calculate zmin and zmax for surface
    vacuum = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin + conf.geometry.box_height * coarse_surf.sizes.zbox);
    bulk = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin - conf.geometry.bulk_height * conf.geometry.latconst);
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

    Surface bulk, coarse_surf, vacuum;
    fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
    if (fail) return 1;

    start_msg(t0, "=== Making big mesh...");
    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality, a - element volume,
    // F - suppress output of faces and edges, B - suppress output of boundary info
    string command = "rnQFBq" + conf.geometry.mesh_quality;
    if (conf.geometry.element_volume != "") command += "a" + conf.geometry.element_volume;
    int err_code = fem_mesh.generate(bulk, coarse_surf, vacuum, command);
    check_return(err_code, "Triangulation failed with error code " + to_string(err_code));
    end_msg(t0);

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = fem_mesh.mark_mesh();
    check_return(fail, "Mesh marking failed!");
    end_msg(t0);

    start_msg(t0, "=== Generating surface faces...");
    err_code = fem_mesh.generate_surface(reader.sizes, "rnQB");
    end_msg(t0);
    check_return(err_code, "Generation of surface faces failed with error code " + to_string(err_code));

    if (conf.smoothing.algorithm != "none" && conf.smoothing.n_steps > 0) {
        start_msg(t0, "=== Smoothing surface faces...");
        fem_mesh.smoothen(conf.smoothing.n_steps, conf.smoothing.lambda_mesh, conf.smoothing.mu_mesh, conf.smoothing.algorithm);
        end_msg(t0);
    }

    fem_mesh.nodes.write("out/tetmesh_nodes.vtk");
    fem_mesh.faces.write("out/trimesh.vtk");
    fem_mesh.elems.write("out/tetmesh.vtk");

    if (conf.run.surface_cleaner) {
        start_msg(t0, "=== Cleaning surface atoms...");
        dense_surf.clean_by_triangles(atom2face, vacuum_surface_interpolator, conf.geometry.latconst);
        end_msg(t0);
        dense_surf.write("out/surface_dense_clean.xyz");
    }

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    fem_mesh.generate_hexahedra();
    end_msg(t0);

    fem_mesh.nodes.write("out/hexmesh_nodes.vtk");
    fem_mesh.quads.write("out/quadmesh.vtk");
    fem_mesh.hexahedra.write("out/hexmesh.vtk");
    fem_mesh.write_separate("out/hexmesh_bulk" + timestep_string + ".vtk", TYPES.BULK);
    stringstream ss; ss << fem_mesh;
    write_verbose_msg(ss.str());

    return 0;
}

// Solve Laplace equation
int Femocs::solve_laplace(const double E0) {
    conf.laplace.E0 = E0;       // reset long-range electric field

    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(E0, conf.tolerance.field_min, conf.tolerance.field_max, conf.geometry.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !laplace_solver.import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_vacuum());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Initializing Laplace solver...");
    laplace_solver.set_applied_efield(-E0);
    laplace_solver.setup_system();
    laplace_solver.assemble_system();
    end_msg(t0);

    stringstream ss; ss << laplace_solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Laplace solver...");
    laplace_solver.solve(conf.laplace.n_phi, conf.laplace.phi_error, true, conf.laplace.ssor_param);
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    fail = vacuum_interpolator->extract_solution(&laplace_solver);
    fail |= vacuum_surface_interpolator->extract_solution(&laplace_solver);
    end_msg(t0);

    check_return(fields.check_limits(vacuum_interpolator->get_solutions()), "Field enhancement is out of limits!");

    vacuum_interpolator->write("out/result_E_phi_vacuum.xyz");
    vacuum_interpolator->write("out/result_E_phi_vacuum.vtk");
    vacuum_surface_interpolator->write("out/result_E_phi_surface.xyz");
    vacuum_surface_interpolator->write("out/result_E_phi_surface.vtk");

    start_msg(t0, "=== Extracting general E and phi...");
    fail = general_interpolator.extract_solution(&laplace_solver);
    end_msg(t0);

    general_interpolator.nodes.write("out/result_E_phi.xyz");
    general_interpolator.lintets.write("out/result_E_phi_linear.vtk");
    general_interpolator.quadtets.write("out/result_E_phi_quad.vtk");

    return fail;
}

// Pick a method to solve heat & continuity equations
int Femocs::solve_heat(const double T_ambient) {
    const double delta_t_MD = 4.05e-15; // in seconds

    if (conf.heating.mode == "stationary") {
        return solve_stationary_heat();
    }
    else if(conf.heating.mode == "transient") {
        double delta_time = delta_t_MD * (timestep - last_full_timestep); //in sec
        int success = solve_transient_heat(delta_time);
        write_verbose_msg("Current and heat advanced for " + to_string(delta_time));
        return success;
    }
    else if (conf.heating.mode == "converge") {
        return solve_converge_heat();
    }

    return 0;
}

// Solve steady-state heat and continuity equations
int Femocs::solve_stationary_heat() {
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
    FieldReader fr(NULL, vacuum_surface_interpolator);
    fr.transfer_elfield(ch_solver, conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
    end_msg(t0);

    fr.write("out/surface_field.xyz");

    start_msg(t0, "=== Running J & T solver...\n");
    ch_solver->set_ambient_temperature(conf.heating.t_ambient);
    double t_error = ch_solver->run_specific(conf.heating.t_error, conf.heating.n_newton, false, "", MODES.VERBOSE, 2, 400, true);
    end_msg(t0);

    ch_solver->output_results("out/result_J_T.vtk");

    check_return(t_error > conf.heating.t_error, "Temperature didn't converge, err=" + to_string(t_error));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator->extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator->write("out/result_J_T.xyz");

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
int Femocs::solve_transient_heat(const double delta_time) {
    static bool first_call = true;
    double multiplier = 1.;

    start_msg(t0, "=== Importing mesh to transient J & T solver...");
    fail = !ch_transient_solver.import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_bulk());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Transfering elfield to J & T solver...");
    FieldReader field_reader(NULL, vacuum_surface_interpolator);
    field_reader.transfer_elfield(ch_transient_solver, conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
    end_msg(t0);
    field_reader.write("out/surface_field.xyz");

    start_msg(t0, "=== Interpolating J & T on face centroids...");

    HeatReader heat_reader(NULL, bulk_surface_interpolator);
    heat_reader.interpolate_2d(ch_transient_solver, conf.heating.t_ambient, false);
    end_msg(t0);
    heat_reader.write("out/surface_temperature.xyz");

    start_msg(t0, "=== Calculating field emission...");
    EmissionReader emission(field_reader, heat_reader, fem_mesh.faces, vacuum_interpolator);
    emission.transfer_emission(ch_transient_solver,
            conf.heating.work_function, conf.heating.Vappl);
    end_msg(t0);
    emission.write("out/surface_emission.xyz");

    if (first_call) {
        start_msg(t0, "=== Setup transient J & T solver...");
        ch_transient_solver.setup_current_system();
        ch_transient_solver.setup_heating_system();
        end_msg(t0);
    }

    start_msg(t0, "=== Calculating current density...");
    ch_transient_solver.assemble_current_system(); // assemble matrix for current density equation; current == electric current
    unsigned int ccg = ch_transient_solver.solve_current();  // ccg == number of current calculation (CG) iterations
    end_msg(t0);

    start_msg(t0, "=== Calculating temperature distribution...\n");
    ch_transient_solver.set_timestep(delta_time);

    ch_transient_solver.assemble_heating_system_euler_implicit();

    unsigned int hcg = ch_transient_solver.solve_heat(); // hcg == number of temperature calculation (CG) iterations
    end_msg(t0);

    ch_transient_solver.output_results_current("out/result_J.vtk");
    ch_transient_solver.output_results_heating("out/result_T.vtk");

    start_msg(t0, "=== Extracting J & T...");
    bulk_surface_interpolator->extract_solution(&ch_transient_solver);
    end_msg(t0);
    bulk_surface_interpolator->write("out/result_J_T.movie");

    first_call = false;
    return 0;
}

// Solve transient heat and continuity until convergence is achieved
int Femocs::solve_converge_heat() {
    static bool first_call = true;

    start_msg(t0, "=== Importing mesh to transient J & T solver...");
    fail = !ch_transient_solver.import_mesh_directly(fem_mesh.nodes.export_dealii(),
            fem_mesh.hexahedra.export_bulk());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Transfering elfield to J & T solver...");
    FieldReader field_reader(NULL, vacuum_surface_interpolator);
    field_reader.transfer_elfield(ch_transient_solver, conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
    end_msg(t0);
    field_reader.write("out/surface_field.xyz");

    if (first_call) {
        start_msg(t0, "=== Setup transient J & T solver...");
        ch_transient_solver.setup_current_system();
        ch_transient_solver.setup_heating_system();
        end_msg(t0);
        first_call = false;
    }

    double current_time = 0.;
    double delta_time = 1.e-12; //in seconds!!
    double multiplier = 1.;

    for (int i = 0; i < 1000; ++i){

        start_msg(t0, "=== Interpolating J & T on face centroids...");
        HeatReader heat_reader(NULL, bulk_surface_interpolator);
        heat_reader.interpolate_2d(ch_transient_solver, conf.heating.t_ambient, false);
        end_msg(t0);
        if (MODES.VERBOSE) heat_reader.write("out/surface_temperature.xyz");

        start_msg(t0, "=== Calculating field emission...");
        EmissionReader emission(field_reader, heat_reader, fem_mesh.faces, vacuum_interpolator);
        emission.set_multiplier(multiplier);
        emission.transfer_emission(ch_transient_solver, conf.heating.work_function,
                conf.heating.Vappl, conf.heating.blunt);
        multiplier = emission.get_multiplier();
        end_msg(t0);
        if (MODES.VERBOSE) emission.write("out/surface_emission.xyz");

        start_msg(t0, "=== Calculating current density...");
        ch_transient_solver.assemble_current_system(); // assemble matrix for electric current density equation
        unsigned int ccg = ch_transient_solver.solve_current();  // ccg == number of current calculation (CG) iterations
        end_msg(t0);

        start_msg(t0, "=== Calculating temperature distribution...\n");
        ch_transient_solver.set_timestep(delta_time);
        ch_transient_solver.assemble_heating_system_euler_implicit();
        unsigned int hcg = ch_transient_solver.solve_heat(); // hcg == number of temperature calculation (CG) iterations
        end_msg(t0);

        start_msg(t0, "=== Extracting J & T...");
        bulk_surface_interpolator->extract_solution(&ch_transient_solver);
        end_msg(t0);
        bulk_surface_interpolator->write("out/result_J_T.movie");

        current_time += delta_time;
        if (MODES.VERBOSE) {
            double max_T = ch_transient_solver.get_max_temperature();
            printf("  i=%d ; dt= %5.3e ps ; t= %5.3e ps ; ccg= %2d ; hcg= %2d ; Tmax= %6.2f \n",
                    i, delta_time * 1.e12, current_time * 1.e12, ccg, hcg, max_T);
        }

        if (max(hcg, ccg) < 120 || hcg < 30)
            delta_time *= 1.25;
        else if (max(hcg, ccg) > 150)
            delta_time /= 1.25;

        if (max(hcg, ccg) < 10) return 0;
    }
    write_silent_msg("WARNING: Heat equation did not converged after 1000 steps.");

    return 0;
}

// Generate artificial nanotip
int Femocs::generate_nanotip(const double height, const double radius, const double resolution) {
    clear_log();

    double res = conf.geometry.latconst;
    if (resolution > 0) res = resolution;

    double r = conf.geometry.radius - res;
    if (radius > 0) r = radius;

    start_msg(t0, "=== Generating nanotip...");
    reader.generate_nanotip(height, r, res);
    reader.calc_coordinations(conf.geometry.nnn);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    reader.write("out/atomreader.ckx");
    return 0;
}

// import atoms from a file
int Femocs::import_atoms(const string& file_name, const int add_noise) {
    clear_log();
    string file_type, fname;

    if (file_name == "") fname = conf.path.infile;
    else fname = file_name;

    file_type = get_file_type(fname);
    require(file_type == "ckx" || file_type == "xyz", "Unknown file type: " + file_type);

    start_msg(t0, "=== Importing atoms...");
    reader.import_file(fname, add_noise);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.tolerance.distance) < conf.tolerance.distance;
    end_msg(t0);

    if (!skip_calculations) {
        if (file_type == "xyz") {
            start_msg(t0, "=== Performing coordination analysis...");
            if (!conf.run.rdf) reader.calc_coordinations(conf.geometry.nnn, conf.geometry.coordination_cutoff);
            else reader.calc_coordinations(conf.geometry.nnn, conf.geometry.latconst, conf.geometry.coordination_cutoff);
            end_msg(t0);

            if (conf.run.rdf) {
                stringstream stream;
                stream << fixed << setprecision(3)
                        << "nnn: " << conf.geometry.nnn << ", latconst: " << conf.geometry.latconst
                        << ", coord_cutoff: " << conf.geometry.coordination_cutoff
                        << ", cluster_cutoff: " << conf.geometry.cluster_cutoff;
                write_verbose_msg(stream.str());
            }

            if (conf.run.cluster_anal) {
                start_msg(t0, "=== Performing cluster analysis...");
                reader.calc_clusters(conf.geometry.nnn, conf.geometry.cluster_cutoff, conf.geometry.coordination_cutoff);
                end_msg(t0);
                reader.check_clusters(1);
            }

            start_msg(t0, "=== Extracting atom types...");
            reader.extract_types(conf.geometry.nnn, conf.geometry.coordination_cutoff);
            end_msg(t0);

        } else {
            start_msg(t0, "=== Calculating coords from atom types...");
            reader.calc_coordinations(conf.geometry.nnn);
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
    skip_calculations = reader.calc_rms_distance(conf.tolerance.distance) < conf.tolerance.distance;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Performing coordination analysis...");
        if (!conf.run.rdf) reader.calc_coordinations(conf.geometry.nnn, conf.geometry.coordination_cutoff, nborlist);
        else reader.calc_coordinations(conf.geometry.nnn, conf.geometry.coordination_cutoff, conf.geometry.latconst, nborlist);
        end_msg(t0);

        if (conf.run.rdf) {
            stringstream stream;
            stream << fixed << setprecision(3)
                    << "nnn: " << conf.geometry.nnn << ", latconst: " << conf.geometry.latconst
                    << ", coord_cutoff: " << conf.geometry.coordination_cutoff
                    << ", cluster_cutoff: " << conf.geometry.cluster_cutoff;
            write_verbose_msg(stream.str());
        }

        if (conf.run.cluster_anal) {
            start_msg(t0, "=== Performing cluster analysis...");
            reader.calc_clusters(conf.geometry.nnn, conf.geometry.cluster_cutoff, conf.geometry.coordination_cutoff, nborlist);
            end_msg(t0);
            reader.check_clusters(1);
        }

        start_msg(t0, "=== Extracting atom types...");
        reader.extract_types(conf.geometry.nnn, conf.geometry.latconst);
        end_msg(t0);
    }

    reader.write("out/atomreader.ckx");
    return 0;
}

// import coordinates and types of atoms
int Femocs::import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    clear_log();
    conf.run.surface_cleaner = false; // disable the surface cleaner for atoms with known types

    start_msg(t0, "=== Importing atoms...");
    reader.import_helmod(n_atoms, x, y, z, types);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

    start_msg(t0, "=== Comparing with previous run...");
    skip_calculations = reader.calc_rms_distance(conf.tolerance.distance) < conf.tolerance.distance;
    end_msg(t0);

    if (!skip_calculations) {
        start_msg(t0, "=== Calculating coordinations from atom types...");
        reader.calc_coordinations(conf.geometry.nnn);
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
    check_return(vacuum_interpolator->size() == 0, "No field to export!");

    fail = false;

    if (skip_calculations)
        write_silent_msg("Using previous electric field!");
    else {
        start_msg(t0, "=== Interpolating E and phi...");
        fields.interpolate_2d(atom2face, dense_surf, true);
//        fields.interpolate(dense_surf, 0, true);
        fail = fields.clean(conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
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
    if (n_atoms < 0 || conf.heating.mode == "none") return 0;
    check_return(bulk_interpolator->size() == 0, "No temperature to export!");

    if (skip_calculations)
        write_silent_msg("Using previous temperature!");
    else {
        start_msg(t0, "=== Interpolating J & T...");
        temperatures.interpolate(reader, conf.heating.t_ambient, false);
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
        const double tot_charge = conf.laplace.E0 * reader.sizes.xbox * reader.sizes.ybox;

        ChargeReader face_charges(vacuum_interpolator); // charges on surface triangles
        face_charges.set_check_params(tot_charge, conf.tolerance.charge_min, conf.tolerance.charge_max);

        start_msg(t0, "=== Calculating face charges...");
        face_charges.calc_charges(fem_mesh, conf.laplace.E0);
        end_msg(t0);
        face_charges.write("out/face_charges.xyz");
        check_return(face_charges.check_limits(), "Face charges are not conserved!");

        face_charges.clean(dense_surf.sizes, conf.geometry.latconst);

        start_msg(t0, "=== Distributing face charges...");
        forces.distribute_charges(fields, face_charges, conf.run.hist_cleaner*conf.geometry.coordination_cutoff,
                conf.smoothing.beta_charge);
        end_msg(t0);

        start_msg(t0, "=== Calculating Voronoi charges & forces...");
        VoronoiMesh voro_mesh;
        int err_code;
        err_code = forces.calc_voronoi_charges(voro_mesh, atom2face, fields, conf.geometry.radius, conf.geometry.latconst, "10.0");
        check_return(err_code, "Generation of Voronoi cells failed with error code " + to_string(err_code));
        end_msg(t0);

        voro_mesh.nodes.write("out/voro_nodes.vtk");
        voro_mesh.voros.write("out/voro_cells.vtk");
        voro_mesh.vfaces.write("out/voro_faces.vtk");
        forces.write("out/forces.movie");
        check_return(face_charges.check_limits(forces.get_interpolations()), "Voronoi charges are not conserved!");
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
    check_return(vacuum_surface_interpolator->size() == 0, "No solution to export!");

    FieldReader fr(NULL, vacuum_surface_interpolator);
    start_msg(t0, "=== Interpolating & exporting surface elfield...");
    fr.interpolate_2d(n_points, x, y, z, false);
    fail = fr.clean(conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
    fr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    fr.write("out/interpolation_E_surf.movie");
    return fields.check_limits(fr.get_interpolations()) || fail;
}

// linearly interpolate electric field at given points
int Femocs::interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator->size() == 0, "No electric field to export!");

    FieldReader fr(vacuum_interpolator);
    start_msg(t0, "=== Interpolating & exporting elfield...");
    fr.interpolate(n_points, x, y, z, false);
    fail = fr.clean(conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
    fr.export_elfield(n_points, Ex, Ey, Ez, Enorm, flag);
    end_msg(t0);

    fr.write("out/interpolation_E.movie");
    return fields.check_limits(fr.get_interpolations()) || fail;
}

// linearly interpolate electric potential at given points
int Femocs::interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
        double* phi, int* flag) {

    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator->size() == 0, "No electric potential to export!");

    FieldReader fr(vacuum_interpolator);
    fr.interpolate(n_points, x, y, z, false);
    fail = fr.clean(conf.geometry.coordination_cutoff, conf.run.hist_cleaner);
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
