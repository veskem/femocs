/*
 * ProjectRunaway.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include <omp.h>

#include "ProjectRunaway.h"
#include "Macros.h"
#include "Tethex.h"
#include "VoronoiMesh.h"


using namespace std;
namespace femocs {

ProjectRunaway::ProjectRunaway(AtomReader &reader, Config &config) :
        GeneralProject(reader, config),
        fail(false), t0(0), timestep(-1), last_full_timestep(0),

        vacuum_interpolator("elfield", "potential"),
        bulk_interpolator("rho", "temperature"),

        surface_fields(&vacuum_interpolator),
        surface_temperatures(&bulk_interpolator),

        phys_quantities(config.heating),
        poisson_solver(NULL, &config.field, &vacuum_interpolator.linhex),
        ch_solver(&phys_quantities, &config.heating),

        emission(&surface_fields, &surface_temperatures, &poisson_solver, &vacuum_interpolator),
        pic_solver(&poisson_solver, &ch_solver, &emission, conf.behaviour.rnd_seed)
{
    fields.set_interpolator(&vacuum_interpolator);
    temperatures.set_interpolator(&bulk_interpolator);
    forces.set_interpolator(&vacuum_interpolator);
    poisson_solver.set_particles(pic_solver.get_particles());

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    phys_quantities.initialize_with_hc_data();
    end_msg(t0);
}

int ProjectRunaway::reinit(int tstep) {
    static bool skip_meshing = true;  // remember the value of skip_meshing
    if (tstep >= 0)
        timestep = tstep;
    else
        ++timestep;

    if (!skip_meshing && MODES.WRITEFILE)
        MODES.WRITEFILE = false;

    if ((conf.behaviour.n_writefile > 0) && (timestep % conf.behaviour.n_writefile == 0))
        MODES.WRITEFILE = true;

    atom2face.clear();

    timestep_string = d2s(timestep);
    timestep_string = "_" + string( max(0.0, 6.0 - timestep_string.length()), '0' ) + timestep_string;

    mesh_changed = false;
    skip_meshing = reader.rms_distance < conf.tolerance.distance;

    write_silent_msg("Running at timestep " + d2s(timestep));
    return skip_meshing;
}

int ProjectRunaway::finalize(double tstart) {
    reader.save_current_run_points(conf.tolerance.distance);

    // In case no transient simulations (time has stayed the same) advance the time
    if (conf.pic.mode == "none" && conf.heating.mode == "none")
        GLOBALS.TIME += conf.behaviour.timestep_fs;

    last_full_timestep = timestep;

    write_silent_msg("Total execution time " + d2s(omp_get_wtime()-tstart, 3));
    return 0;
}

int ProjectRunaway::run(const int timestep) {
    return run(conf.field.E0, timestep);
}

int ProjectRunaway::run(const double elfield, const int tstep) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (reinit(tstep))
        write_verbose_msg("Atoms haven't moved significantly, "
                + d2s(reader.rms_distance, 3) + " < " + d2s(conf.tolerance.distance, 3)
                + "! Previous mesh will be used!");

    else if (generate_mesh())
        return process_failed("Mesh generation failed!");

    check_return(GLOBALS.TIME == 0 && !mesh_changed, "First meshing failed! Terminating...");

    if (mesh_changed && prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    //***** Run FEM solvers *****

    if (run_field_solver())
        return process_failed("Running field solver in a " + conf.field.solver + " mode failed!");

    if (run_heat_solver())
        return process_failed("Running heat solver in a " + conf.heating.mode + " mode failed!");

    //***** Prepare for data export and next run *****

    if (prepare_export())
        return process_failed("Interpolating solution on atoms failed!");

    finalize(tstart);
    return 0;
}

int ProjectRunaway::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum) {
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
    double box_height = max(conf.geometry.latconst, coarse_surf.sizes.zbox) * conf.geometry.box_height;
    vacuum = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin + box_height);
    bulk = Surface(coarse_surf.sizes, coarse_surf.sizes.zmin - conf.geometry.bulk_height * conf.geometry.latconst);
    reader.resize_box(coarse_surf.sizes.xmin, coarse_surf.sizes.xmax,
            coarse_surf.sizes.ymin, coarse_surf.sizes.ymax,
            bulk.sizes.zmin, vacuum.sizes.zmax);
    end_msg(t0);

    bulk.write("out/bulk.xyz");
    vacuum.write("out/vacuum.xyz");

    return 0;
}

int ProjectRunaway::generate_mesh() {
    if (conf.path.mesh_file != "")
        return read_mesh();

    new_mesh->clear();

    Surface bulk, coarse_surf, vacuum;
    fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
    if (fail) return 1;

    start_msg(t0, "=== Making big mesh...");
    // r - reconstruct, n(n) - output tet neighbour list (and tri-tet connection),
    // Q - quiet, q - mesh quality, a - element volume, E - suppress output of elements
    // F - suppress output of faces and edges, B - suppress output of boundary info
    string command = "rQFBq" + conf.geometry.mesh_quality;
    if (conf.geometry.element_volume != "") command += "a" + conf.geometry.element_volume;
    int err_code = new_mesh->generate(bulk, coarse_surf, vacuum, command);
    check_return(err_code, "Triangulation failed with error code " + d2s(err_code));
    end_msg(t0);

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = new_mesh->mark_mesh();
    check_return(fail, "Mesh marking failed!");
    end_msg(t0);

    start_msg(t0, "=== Generating surface faces...");
    err_code = new_mesh->generate_surface(reader.sizes, "rQB", "rQnn");
    end_msg(t0);
    check_return(err_code, "Generation of surface faces failed with error code " + d2s(err_code));

    if (conf.smoothing.algorithm != "none" && conf.smoothing.n_steps > 0) {
        start_msg(t0, "=== Smoothing surface faces...");
        new_mesh->smoothen(conf.smoothing.n_steps, conf.smoothing.lambda_mesh, conf.smoothing.mu_mesh, conf.smoothing.algorithm);
        end_msg(t0);
    }

    // has to be separate to ensure that data is for sure calculated
    new_mesh->tris.calc_appendices();

    new_mesh->nodes.write("out/tetmesh_nodes.vtk");
    new_mesh->tris.write("out/trimesh.vtk");
    new_mesh->tets.write("out/tetmesh.vtk");

    if (conf.run.surface_cleaner) {
        start_msg(t0, "=== Cleaning surface atoms...");
        dense_surf.clean_by_triangles(atom2face, vacuum_interpolator, new_mesh, conf.geometry.latconst);
        end_msg(t0);
        dense_surf.write("out/surface_dense_clean.xyz");
    }

    start_msg(t0, "=== Converting tetrahedra to hexahedra...");
    new_mesh->generate_hexahedra();
    end_msg(t0);

    new_mesh->nodes.write("out/hexmesh_nodes.vtk");
    new_mesh->quads.write("out/quadmesh.vtk");
    new_mesh->hexs.write("out/hexmesh.vtk");
    new_mesh->write_separate("out/hexmesh_bulk" + timestep_string + ".vtk", TYPES.BULK);

    // update mesh pointers
    mesh = new_mesh;
    mesh_changed = true;

    static bool odd_run = true;
    if (odd_run) new_mesh = &mesh2;
    else new_mesh = &mesh1;
    odd_run = !odd_run;

    write_verbose_msg(mesh->to_str());
    return 0;
}

int ProjectRunaway::read_mesh() {
    start_msg(t0, "=== Reading mesh from file...");
    int err_code = new_mesh->read(conf.path.mesh_file, "rQnn");
    end_msg(t0);
    check_return(err_code, "Reading mesh failed with error code " + d2s(err_code));

    new_mesh->nodes.write("out/hexmesh_nodes.vtk");
    new_mesh->nodes.write("out/hexmesh_nodes.xyz");
    new_mesh->tris.write("out/trimesh.vtk");
    new_mesh->quads.write("out/quadmesh.vtk");
    new_mesh->tets.write("out/tetmesh.vtk");
    new_mesh->hexs.write("out/hexmesh.vtk");
    new_mesh->write_separate("out/hexmesh_bulk" + timestep_string + ".vtk", TYPES.BULK);

    mesh = new_mesh;
    mesh_changed = true;

    write_verbose_msg(mesh->to_str());
    return 0;
}

int ProjectRunaway::prepare_solvers() {
    start_msg(t0, "=== Importing mesh to Poisson solver...");
    fail = !poisson_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);

    vacuum_interpolator.initialize(mesh, poisson_solver);
    vacuum_interpolator.lintet.narrow_search_to(TYPES.VACUUM);

    if (conf.field.solver == "poisson" || conf.heating.mode != "none") {
        start_msg(t0, "=== Importing mesh to J & T solver...");
        fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_bulk());
        check_return(fail, "Importing bulk mesh to Deal.II failed!");
        end_msg(t0);

        bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient);
        bulk_interpolator.lintet.narrow_search_to(TYPES.BULK);

        surface_fields.set_preferences(false, 2, 3);
        surface_temperatures.set_preferences(false, 2, 3);
    }

    return 0;
}

int ProjectRunaway::prepare_export() {
    if (dense_surf.size() == 0) {
        write_silent_msg("No atoms for solution interpolation!");
        return 0;
    }

    start_msg(t0, "=== Interpolating E and phi...");
    fields.set_preferences(true, 2, 1);
    fields.interpolate(dense_surf);
    end_msg(t0);

    fields.write("out/fields.movie");
    check_return(fields.check_limits(), "Field enhancement is out of limits!");

    if (conf.heating.mode != "none") {
        start_msg(t0, "=== Interpolating J & T...");
        temperatures.set_preferences(true, 3, conf.behaviour.interpolation_rank);
        temperatures.interpolate(reader);
        end_msg(t0);

        // TODO implement reasonable temperature limit check
        temperatures.write("out/temperatures.movie");
    }

    if (conf.force.mode != "none") {
        // analytical total charge without epsilon0 (will be added in ChargeReader)
        const double tot_charge = conf.field.E0 * reader.sizes.xbox * reader.sizes.ybox;

        ChargeReader face_charges(&vacuum_interpolator); // charges on surface triangles
        face_charges.set_check_params(tot_charge, conf.tolerance.charge_min, conf.tolerance.charge_max);

        start_msg(t0, "=== Calculating face charges...");
        face_charges.calc_charges(*mesh, conf.field.E0);
        end_msg(t0);

        face_charges.write("out/face_charges.xyz");
        check_return(face_charges.check_limits(), "Face charges are not conserved!");

        start_msg(t0, "=== Distributing face charges...");
        // Remove the atoms and their solutions outside the box
        face_charges.clean(dense_surf.sizes, conf.geometry.latconst);
        forces.distribute_charges(fields, face_charges, conf.smoothing.beta_charge);
        end_msg(t0);

        start_msg(t0, "=== Generating Voronoi cells...");
        VoronoiMesh voro_mesh;
        int err_code = forces.calc_voronois(voro_mesh, atom2face, conf.geometry.radius, conf.geometry.latconst, "10.0");
        end_msg(t0);

        check_return(err_code, "Generation of Voronoi cells failed with error code " + d2s(err_code));
        voro_mesh.nodes.write("out/voro_nodes.vtk");
        voro_mesh.voros.write("out/voro_cells.vtk");
        voro_mesh.vfaces.write("out/voro_faces.vtk");

        if (conf.force.mode == "lorentz") {
            start_msg(t0, "=== Calculating Lorentz force...");
            forces.calc_charge_and_lorentz(voro_mesh, fields);
        } else {
            start_msg(t0, "=== Calculating Lorentz & Coulomb force...");
            forces.calc_charge_and_lorentz(voro_mesh, fields);
            forces.calc_coulomb(conf.geometry.charge_cutoff);
        }
        end_msg(t0);

        forces.write("out/forces.movie");
        check_return(face_charges.check_limits(forces.get_interpolations()), "Voronoi charges are not conserved!");
    }

    return 0;
}

int ProjectRunaway::run_field_solver() {
    if (conf.field.solver == "poisson") {
        if (conf.pic.mode == "transient")
            return solve_pic(conf.behaviour.timestep_fs, mesh_changed);
        else if (conf.pic.mode == "converge")
            return converge_pic(1.e4);
        else
            check_return(false, "Invalid PIC mode: " + conf.pic.mode);
    }

    if (mesh_changed && (conf.field.solver == "laplace" || conf.pic.mode == "none"))
        return solve_laplace(conf.field.E0, conf.field.V0);

    return 0;
}

int ProjectRunaway::run_heat_solver() {
    int ccg, hcg;

    if (conf.heating.mode == "converge")
        return converge_heat(conf.heating.t_ambient);

    if (mesh_changed && conf.heating.mode == "transient")
        return solve_heat(conf.heating.t_ambient, GLOBALS.TIME - last_heat_time, true, ccg, hcg);

    return 0;
}

int ProjectRunaway::solve_laplace(double E0, double V0) {
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(conf, conf.geometry.radius, mesh->tris.stat.zbox, mesh->nodes.stat.zbox);

    start_msg(t0, "=== Initializing Laplace solver...");
    poisson_solver.setup(-E0, V0);
    poisson_solver.assemble_laplace(true);
    end_msg(t0);

    write_verbose_msg(poisson_solver.to_str());

    start_msg(t0, "=== Running Laplace solver...");
    int ncg = poisson_solver.solve();
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    vacuum_interpolator.extract_solution(poisson_solver);
    end_msg(t0);

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintet.write("out/result_E_phi.vtk");

    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()),
            "Field enhancement is out of limits!");
    return 0;
}

int ProjectRunaway::solve_pic(double advance_time, bool reinit) {
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(conf, conf.geometry.radius, mesh->tris.stat.zbox, mesh->nodes.stat.zbox);

    int n_pic_steps = ceil(advance_time / conf.pic.dt_max); // dt_main = delta_t_MD converted to [fs]
    double dt_pic = advance_time / n_pic_steps;

    // vector holding point to cell mapping to make interpolation faster
    vector<int> point2cell;

    if (reinit) {
        start_msg(t0, "=== Initializing Poisson solver...");
        poisson_solver.setup(-conf.field.E0, conf.field.V0);
        end_msg(t0);
        write_verbose_msg(poisson_solver.to_str());
    }

    pic_solver.set_params(conf.pic, dt_pic, mesh->nodes.stat);
    surface_temperatures.interpolate(ch_solver);
    surface_fields.store_points(ch_solver);
    emission.initialize(mesh);

    start_msg(t0, "=== Running PIC for delta time " + d2s(advance_time, 2) + " fs\n");
    for (int i = 0; i < n_pic_steps; ++i) {
        // advance super particles
        int n_lost = pic_solver.update_positions();

        // assemble and solve Poisson equation
        poisson_solver.assemble_poisson(i == 0, write_time());
        int n_cg_steps = poisson_solver.solve();

        // must be after solving Poisson and before updating velocities
        vacuum_interpolator.extract_solution(poisson_solver);

        // update velocities of super particles and collide them
        pic_solver.update_velocities();
        pic_solver.collide_particles();

        // update field on the surface
        surface_fields.calc_interpolation(point2cell);

        // calculate field emission and inject electrons
        emission.calc_emission(conf.emission, conf.field.V0);
        int n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);
        
        write_results();
        if (MODES.VERBOSE)
            printf("  t=%.2e fs, #CG=%d, Fmax=%.3f V/A, Itot=%.3e A, #el inj|del|tot=%d|%d|%d\n",
                    GLOBALS.TIME, n_cg_steps, emission.global_data.Fmax, emission.global_data.I_tot,
                    n_injected, n_lost, pic_solver.get_n_electrons());

        GLOBALS.TIME += dt_pic;
    }
    end_msg(t0);

    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()),
            "Field enhancement is out of limits!");
    return 0;

    // TODO Save ions and neutrals that are inbound on the MD domain somewhere where the MD can find them
}

int ProjectRunaway::converge_pic(double max_time) {
    double time_window; //time window to check convergence
    int i_max; //window iterations
    if (max_time < conf.pic.dt_max * 32) {
        time_window = max_time;
        i_max = 1;
    } else {
        i_max =  ceil(max_time / (25 * conf.pic.dt_max));
        time_window = max_time / i_max;
    }

    double I_mean, I_mean_prev;

    start_msg(t0, "=== Converging PIC with time window " + d2s(time_window, 2) + " fs\n");
    for (int i = 0; i < i_max; ++i) {
        I_mean_prev = emission.get_mean_current();
        solve_pic(time_window, i==0);
        I_mean = emission.get_mean_current();

        double err = (I_mean - I_mean_prev) / I_mean;
        if (MODES.VERBOSE)
            printf("  i=%d, I_mean=%e, error=%e\n", i, I_mean, err);

        if (fabs(err) < 0.05 && fabs(err) < conf.pic.convergence * I_std / I_mean)
            break;
    }
    return 0;
}

int ProjectRunaway::solve_heat(double T_ambient, double delta_time, bool full_run, int& ccg, int& hcg) {
    if (full_run)
        ch_solver.setup(T_ambient);

    // Calculate field emission in case not ready from pic
    if (conf.pic.mode == "none" || conf.field.solver == "laplace") {
        start_msg(t0, "=== Calculating electron emission...");
        if (full_run) surface_fields.interpolate(ch_solver);
        surface_temperatures.interpolate(ch_solver);
        emission.initialize(mesh);
        emission.calc_emission(conf.emission, conf.field.V0);
        end_msg(t0);
    }

    emission.export_emission(ch_solver);

    start_msg(t0, "=== Calculating current density...");
    ch_solver.current.assemble();
    ccg = ch_solver.current.solve();
    end_msg(t0);
    write_verbose_msg("#CG steps: " + d2s(ccg));

    start_msg(t0, "=== Calculating temperature distribution...");
    ch_solver.heat.assemble(delta_time * 1.e-15); // caution!! ch_solver internal time in sec
    hcg = ch_solver.heat.solve();
    end_msg(t0);
    write_verbose_msg("#CG steps: " + d2s(hcg));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    if(conf.pic.mode != "transient")
        GLOBALS.TIME += delta_time;

    write_results();

    last_heat_time = GLOBALS.TIME;
    return 0;
}

int ProjectRunaway::converge_heat(double T_ambient) {
    const int max_steps = 1000;
    double delta_time = conf.heating.delta_time;
    int ccg, hcg, step, error;

//    mesh_changed = false; // run convergence always with the same mesh
    bool global_verbosity = MODES.VERBOSE;

    start_msg(t0, "=== Converging heat...\n");
//    MODES.VERBOSE = false;

    for (step = 0; step < max_steps; ++step) {

        // advance heat and current system for delta_time
        error = solve_heat(conf.heating.t_ambient, delta_time, step == 0, ccg, hcg);
        if (error) return error;

        // modify the advance time depending on how slowly the solution is changing
        if (conf.pic.mode == "none" || conf.pic.mode == "converge")
            GLOBALS.TIME += delta_time;

        if (hcg < (ccg - 10)) // heat changed too little?
            delta_time *= 1.25;
        else if (hcg > (ccg + 10)) // heat changed too much?
            delta_time /= 1.25;

        // write debug data
        if (global_verbosity)
            printf( "t=%.2e ps, dt=%.2e ps, Tmax=%.2fK\n",
                    GLOBALS.TIME * 1.e-3, delta_time * 1.e-3, ch_solver.heat.max_solution() );
        write_results();

        // check if the result has converged
        if (max(hcg, ccg) < 10) break;

        // update field - advance PIC for delta time
        if (conf.pic.mode == "transient")
            error = solve_pic(delta_time, false);
        else if (conf.pic.mode == "converge")
            error = converge_pic(delta_time);
        if (error) return error;

    }

    MODES.VERBOSE = global_verbosity;
    end_msg(t0);

    check_return(step < max_steps, "Failed to converge heat equation after " + d2s(max_steps) + " steps!");
    return 0;
}

int ProjectRunaway::write_results(){

    if (!write_time()) return 1;

    vacuum_interpolator.extract_solution(poisson_solver);
    vacuum_interpolator.nodes.write("out/result_E_phi.movie");
    vacuum_interpolator.linhex.write("out/result_E_phi.vtk");

    if (conf.pic.mode != "none"){
        emission.write("out/emission.dat");
        emission.write("out/emission.movie");
        pic_solver.write("out/electrons.movie");
        surface_fields.write("out/surface_fields.movie");
        vacuum_interpolator.extract_charge_density(poisson_solver);
        vacuum_interpolator.nodes.write("out/result_E_charge.movie");
    }

    if (conf.pic.mode == "none" && emission.atoms.size() > 0){
        emission.write("out/emission.dat");
        emission.write("out/surface_emission.movie");
    }

    if (conf.heating.mode != "none"){
        bulk_interpolator.nodes.write("out/result_J_T.movie");
        bulk_interpolator.lintet.write("out/result_J_T.vtk");
    }

    last_write_time = GLOBALS.TIME;
    return 0;
}

int ProjectRunaway::force_output() {
    if (conf.behaviour.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    mesh->hexs.write("out/hexmesh_err.vtk");
    mesh->tets.write("out/tetmesh_err.vtk");
    mesh->tris.write("out/trimesh_err.vtk");

    vacuum_interpolator.nodes.write("out/result_E_phi_err.xyz");
    vacuum_interpolator.lintet.write("out/result_E_phi_err.vtk");

    if (bulk_interpolator.nodes.size() > 0) {
        if (conf.heating.mode == "transient") {
            bulk_interpolator.nodes.write("out/result_J_T_err.xyz");
            ch_solver.current.write("out/result_J_err.vtk");
            ch_solver.heat.write("out/result_T_err.vtk");
        }
    }

    return 0;
}

int ProjectRunaway::export_results(const int n_points, const string &cmd, double* data) {
    if (n_points <= 0) return 0;

    if (cmd == LABELS.elfield.second || cmd == LABELS.elfield_norm.second || cmd == LABELS.potential.second)
        return fields.export_results(n_points, cmd, false, data);

    if (cmd == LABELS.rho.second || cmd == LABELS.rho_norm.second || cmd == LABELS.temperature.second)
        return temperatures.export_results(n_points, cmd, false, data);

    if (cmd == LABELS.force.second || cmd == LABELS.force_norm.second || cmd == LABELS.charge.second)
        return forces.export_results(n_points, cmd, true, data);

    require(false, "Unimplemented type of export data: " + cmd);
    return 1;
}

int ProjectRunaway::interpolate_results(const int n_points, const string &cmd, const bool on_surface,
        const double* x, const double* y, const double* z, double* data, int* flag) {

    // location of interpolation; 2-on surface, 3-in space
    int dim = 3;
    if (on_surface) dim = 2;

    if (cmd == LABELS.elfield.second || cmd == LABELS.elfield_norm.second || cmd == LABELS.potential.second) {
        fields.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        return fields.interpolate_results(n_points, cmd, x, y, z, data);
    }

    if (cmd == LABELS.rho.second || cmd == LABELS.rho_norm.second || cmd == LABELS.temperature.second) {
        temperatures.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        return temperatures.interpolate_results(n_points, cmd, x, y, z, data);
    }

    if (cmd == LABELS.force.second || cmd == LABELS.force_norm.second || cmd == LABELS.charge.second) {
        forces.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        return forces.interpolate_results(n_points, cmd, x, y, z, data);
    }

    require(false, "Unimplemented type of interpolation data: " + cmd);
    return 1;
}

} /* namespace femocs */
