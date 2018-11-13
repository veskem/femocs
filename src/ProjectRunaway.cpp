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
        fail(false), t0(0), mesh_changed(false),
		last_heat_time(0), last_write_time(0),
		last_write_ts(0), last_mesh_write_ts(0),

        vacuum_interpolator("elfield", "elfield_norm", "potential"),
        bulk_interpolator("rho", "potential", "temperature"),

        fields(&vacuum_interpolator),
        temperatures(&bulk_interpolator),
        forces(&vacuum_interpolator),

        surface_fields(&vacuum_interpolator),
        surface_temperatures(&bulk_interpolator),
        heat_transfer(&bulk_interpolator),

        phys_quantities(config.heating),
        poisson_solver(NULL, &config.field, &vacuum_interpolator.linhex),
        ch_solver(&phys_quantities, &config.heating),

        emission(&surface_fields, &surface_temperatures, &poisson_solver, &vacuum_interpolator),
        pic_solver(&poisson_solver, &ch_solver, &emission, conf.behaviour.rnd_seed)
{
    poisson_solver.set_particles(pic_solver.get_particles());
    temperatures.set_params(config);
    ch_solver.current.set_bcs(emission.get_current_densities());
    ch_solver.heat.set_bcs(emission.get_nottingham());

    // Initialise heating module
    start_msg(t0, "Reading physical quantities");
    phys_quantities.initialize_with_hc_data();
    end_msg(t0);
}

int ProjectRunaway::reinit(int tstep, double time) {
    GLOBALS.TIMESTEP++;

    timestep_string = d2s(GLOBALS.TIMESTEP);
    timestep_string = "_" + string( max(0.0, 6.0 - timestep_string.length()), '0' ) + timestep_string;

    mesh_changed = false;
    bool skip_meshing = reader.get_rmsd() < conf.tolerance.distance;

    MODES.WRITEFILE = conf.behaviour.n_writefile > 0 && (
            last_write_ts == 0 ||
            last_mesh_write_ts == 0 ||
            (GLOBALS.TIMESTEP - last_write_ts) >= conf.behaviour.n_writefile
            );

    write_silent_msg("Running at timestep=" + d2s(GLOBALS.TIMESTEP) + ", time=" + d2s(GLOBALS.TIME, 2) + " fs");
    return skip_meshing;
}

int ProjectRunaway::finalize(double tstart, double time) {
    if (conf.field.mode == "laplace")
        GLOBALS.TIME += conf.behaviour.timestep_fs;
    reader.save_current_run_points(conf.tolerance.distance);

    write_silent_msg("Total execution time " + d2s(omp_get_wtime()-tstart, 3));

    // determine whether log file will be written on next timestep
    if (conf.behaviour.n_write_log > 0)
        MODES.WRITELOG = (GLOBALS.TIMESTEP+1)%conf.behaviour.n_write_log == 0;

    // update file output counters
    if (MODES.WRITEFILE) {
        last_write_ts = GLOBALS.TIMESTEP;
        if (mesh_changed) last_mesh_write_ts = GLOBALS.TIMESTEP;
    }

    return 0;
}

int ProjectRunaway::process_failed(const string &msg) {
    write_verbose_msg(msg);
//    force_output();
    GLOBALS.TIME += conf.behaviour.timestep_fs;
    if (conf.behaviour.n_write_log > 0)
        MODES.WRITELOG = (GLOBALS.TIMESTEP+1)%conf.behaviour.n_write_log == 0;
    return 1;
}

int ProjectRunaway::run(const int timestep, const double time) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (reinit(timestep, time)) {
        write_verbose_msg("Atoms haven't moved significantly, "
                + d2s(reader.get_rmsd(), 3) + " < " + d2s(conf.tolerance.distance, 3)
                + "! Previous mesh will be used!");
        dense_surf.update_positions(reader);
    }

    else if (generate_mesh())
        return process_failed("Mesh generation failed!");

    if (mesh_changed && prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    //***** Run FEM solvers *****

    if (run_field_solver())
        return process_failed("Running field solver in a " + conf.field.mode + " mode failed!");

    if (run_heat_solver())
        return process_failed("Running heat solver in a " + conf.heating.mode + " mode failed!");

    //***** Prepare for data export and next run *****

    // Must be here to ensure old mesh is not deleted too early nor too late
    update_mesh_pointers();

    if (prepare_export())
        return process_failed("Interpolating solution on atoms failed!");

    return finalize(tstart, time);
}

void ProjectRunaway::update_mesh_pointers() {
    static bool odd_run = true;
    if (mesh_changed) {
        if (odd_run) new_mesh = &mesh2;
        else new_mesh = &mesh1;
        odd_run = !odd_run;
    }
}

int ProjectRunaway::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum) {
    start_msg(t0, "Extracting surface");
    reader.extract(dense_surf, TYPES.SURFACE);
    end_msg(t0);
    dense_surf.write("out/surface_dense.xyz");

    if (GLOBALS.TIME == 0) {
        start_msg(t0, "Extending surface");
        dense_surf.extend(extended_surf, conf);
        end_msg(t0);
        extended_surf.write("out/surface_extension.xyz");
    }

    start_msg(t0, "Coarsening surface");
    dense_surf.generate_boundary_nodes(bulk, coarse_surf, vacuum, extended_surf, conf, GLOBALS.TIME==0);
    end_msg(t0);

    if (MODES.VERBOSE)
        printf("  #extended=%d, #coarse=%d, #dense=%d\n", extended_surf.size(), coarse_surf.size(), dense_surf.size());

    coarse_surf.write("out/surface_coarse.xyz");
    bulk.write("out/bulk.xyz");
    vacuum.write("out/vacuum.xyz");

    return 0;
}

int ProjectRunaway::generate_mesh() {
    if (conf.path.mesh_file != "") {
        start_msg(t0, "Reading mesh from file");
        fail = new_mesh->read(conf.path.mesh_file, "rQnn");
    } else {
        Surface bulk, coarse_surf, vacuum;
        fail = generate_boundary_nodes(bulk, coarse_surf, vacuum);
        check_return(fail, "Generation of mesh generator nodes failed!");

        start_msg(t0, "Generating vacuum & bulk mesh");
        fail = new_mesh->generate(bulk, coarse_surf, vacuum, conf);
    }
    end_msg(t0);
    if (fail) return 1;

    if (conf.run.surface_cleaner && dense_surf.size() > 0) {
        start_msg(t0, "Cleaning surface atoms");
        dense_surf.clean_by_triangles(vacuum_interpolator, new_mesh, conf.geometry.latconst);
        end_msg(t0);
    }

    MODES.WRITEFILE |= conf.behaviour.n_writefile > 0 &&
            (GLOBALS.TIMESTEP - last_mesh_write_ts) >= conf.behaviour.n_writefile;

    new_mesh->nodes.write("out/hexmesh_nodes.vtk");
    new_mesh->tris.write("out/trimesh.vtk");
    new_mesh->quads.write("out/quadmesh.vtk");
    new_mesh->tets.write("out/tetmesh.vtk");
    new_mesh->hexs.write("out/hexmesh.vtk");
    new_mesh->write_separate("out/hexmesh_bulk" + timestep_string + ".vtk", TYPES.BULK);

    // update mesh pointers
    mesh = new_mesh;
    mesh_changed = true;

    write_verbose_msg(mesh->to_str());
    return 0;
}

int ProjectRunaway::prepare_solvers() {
    start_msg(t0, "Importing vacuum mesh to Deal.II");
    fail = !poisson_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);

    if (conf.field.mode != "laplace" || conf.heating.mode != "none") {
        start_msg(t0, "Importing bulk mesh to Deal.II");
        fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_bulk());
        check_return(fail, "Importing bulk mesh to Deal.II failed!");
        end_msg(t0);

        // in case of first run, give initial values to the temperature
        // interpolator and set interpolation preferences
        if (bulk_interpolator.nodes.size() == 0) {
            bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient, TYPES.BULK);
            surface_fields.set_preferences(false, 2, 3);
            surface_temperatures.set_preferences(false, 2, 3);
            heat_transfer.set_preferences(false, 3, 1);
        }

        // setup ch_solver here as pic_solver might intialize bulk_interpolator
        start_msg(t0, "Setup current & heat solvers");
        ch_solver.setup(conf.heating.t_ambient);
        heat_transfer.interpolate_dofs(ch_solver);  // Transfer previous temperatures to new mesh
        end_msg(t0);
    }

    return 0;
}

int ProjectRunaway::prepare_export() {
    if (dense_surf.size() == 0) {
        write_silent_msg("No atoms for solution interpolation!");
        return 0;
    }

    start_msg(t0, "Interpolating E & phi");
    if (mesh_changed) {
        fields.set_preferences(true, 2, 1);
        fields.interpolate(dense_surf);
    } else {
        fields.update_positions(dense_surf);
        fields.calc_interpolation();
    }
    end_msg(t0);

    fields.write("out/fields.movie");
    check_return(fields.check_limits(), "Field enhancement is out of limits!");

    if (conf.heating.mode != "none") {
        start_msg(t0, "Interpolating J & T");
        if (mesh_changed) {
            temperatures.set_preferences(false, 3, conf.behaviour.interpolation_rank);
            temperatures.interpolate(reader);
            temperatures.precalc_berendsen_long();
        } else {
            temperatures.update_positions(reader);
            temperatures.calc_interpolation();
        }
        end_msg(t0);

        temperatures.write("out/temperatures.xyz");
        // TODO implement reasonable temperature limit check
    }

    int retval = 0;
    if (conf.force.mode != "none") {
        if (mesh_changed)
            retval = solve_force();
        else {
            start_msg(t0, "Recalculating forces");
            forces.recalc_lorentz(fields);
            if (conf.force.mode == "all")
                forces.calc_coulomb(conf.geometry.charge_cutoff);

            end_msg(t0);
        }

        forces.write("out/forces.movie");
    }

    return retval;
}

int ProjectRunaway::run_field_solver() {
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    if (mesh_changed)
        fields.set_check_params(conf, conf.geometry.radius, mesh->tris.stat.zbox, mesh->nodes.stat.zbox);

    if (conf.field.mode != "laplace")
        return solve_pic(conf.behaviour.timestep_fs, mesh_changed);

    else if (mesh_changed)
        return solve_laplace(conf.field.E0, conf.field.V0);

    return 0;
}

int ProjectRunaway::run_heat_solver() {
    int ccg, hcg;

    double delta_time = GLOBALS.TIME - last_heat_time;
    bool b1 = delta_time >= conf.heating.delta_time;
    if (conf.heating.mode == "transient" && (mesh_changed || b1))
        return solve_heat(conf.heating.t_ambient, max(delta_time, conf.behaviour.timestep_fs), mesh_changed, ccg, hcg);

    return 0;
}

int ProjectRunaway::solve_force() {
    // analytical total charge without epsilon0 (will be added in ChargeReader)
    const double tot_charge = conf.field.E0 * mesh->nodes.stat.xbox * mesh->nodes.stat.ybox;

    ChargeReader face_charges(&vacuum_interpolator); // charges on surface triangles
    face_charges.set_check_params(tot_charge, conf.tolerance.charge_min, conf.tolerance.charge_max);

    start_msg(t0, "Calculating face charges");
    face_charges.calc_charges(*mesh, conf.field.E0);
    end_msg(t0);

    face_charges.write("out/face_charges.xyz");
    check_return(face_charges.check_limits(), "Face charges are not conserved!");

    start_msg(t0, "Distributing face charges");
    // Remove the atoms and their solutions outside the box
    face_charges.clean(dense_surf.sizes, conf.geometry.latconst);
    forces.distribute_charges(fields, face_charges, conf.smoothing.beta_charge);
    end_msg(t0);

    start_msg(t0, "Generating Voronoi cells");
    VoronoiMesh voro_mesh;
    int err_code = forces.calc_voronois(voro_mesh, conf.geometry, "10.0");
    end_msg(t0);

    check_return(err_code, "Generation of Voronoi cells failed with error code " + d2s(err_code));
    voro_mesh.nodes.write("out/voro_nodes.vtk");
    voro_mesh.voros.write("out/voro_cells.vtk");
    voro_mesh.vfaces.write("out/voro_faces.vtk");

    if (conf.force.mode == "lorentz") {
        start_msg(t0, "Calculating Lorentz force");
        forces.calc_charge_and_lorentz(voro_mesh, fields);
    } else {
        start_msg(t0, "Calculating Lorentz & Coulomb force");
        forces.calc_charge_and_lorentz(voro_mesh, fields);
        forces.calc_coulomb(conf.geometry.charge_cutoff);
    }
    end_msg(t0);

    forces.write("out/forces.movie");
    check_return(face_charges.check_limits(forces.get_interpolations()), "Voronoi charges are not conserved!");

    return 0;
}

int ProjectRunaway::solve_laplace(double E0, double V0) {
    start_msg(t0, "Initializing Laplace solver");
    poisson_solver.setup(-E0, V0);
    poisson_solver.assemble(true, false);
    end_msg(t0);

    write_verbose_msg(poisson_solver.to_str());

    start_msg(t0, "Running Laplace solver");
    int ncg = poisson_solver.solve();
    end_msg(t0);
    check_return(ncg < 0, "Field solver did not complete normally,"
            " #CG=" + d2s(abs(ncg)) + "/" + d2s(conf.field.n_cg));

    start_msg(t0, "Extracting E & phi");
    vacuum_interpolator.initialize(mesh, poisson_solver, 0, TYPES.VACUUM);
    vacuum_interpolator.extract_solution(poisson_solver, conf.run.field_smoother);
    end_msg(t0);

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintet.write("out/result_E_phi.vtk");

    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()),
            "Field enhancement is out of limits!");
    return 0;
}

int ProjectRunaway::solve_pic(double advance_time, bool full_run) {
    int n_pic_steps = ceil(advance_time / conf.pic.dt_max);
    double dt_pic = advance_time / n_pic_steps;
    pic_solver.set_params(conf.pic, dt_pic, mesh->nodes.stat);

    if (full_run) {
        start_msg(t0, "Initializing Poisson solver");
        poisson_solver.setup(-conf.field.E0, conf.field.V0);
        end_msg(t0);
        write_verbose_msg(poisson_solver.to_str());
    }

    initalize_pic_emission(full_run);

    start_msg(t0, "=== Running PIC for delta time = "
            + d2s(n_pic_steps) + "*" + d2s(dt_pic)+" = " + d2s(advance_time, 2) + " fs\n");
    int n_lost, n_cg, n_injected;
    for (int i = 0; i < n_pic_steps; ++i) {
        make_pic_step(n_lost, n_cg, n_injected, full_run && i==0, write_time());

        if (n_injected > 50000) {
            write_verbose_msg("WARNING: too many injected SP-s, " + d2s(n_injected) + ". Check the SP weight.");
            return 1;
        }

        if (MODES.VERBOSE) {
            printf("  t=%.2e fs, #CG=%d, Fmax=%.3f V/A, Itot=%.3e A, #el inj|del|tot=%d|%d|%d\n",
                    GLOBALS.TIME, n_cg, emission.global_data.Fmax, emission.global_data.I_tot,
                    n_injected, n_lost, pic_solver.get_n_electrons());
        }

        GLOBALS.TIME += dt_pic;
    }
    end_msg(t0);

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintet.write("out/result_E_phi.vtk");

//    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()),
//            "Field enhancement is out of limits!");
    return 0;
}

void ProjectRunaway::initalize_pic_emission(bool full_run) {
    static bool make_intermediate_step = false;

    start_msg(t0, "Calculating surface temperature");
    if (full_run) {
        ch_solver.export_surface_centroids(surface_fields);
        surface_temperatures.interpolate(ch_solver);

        vacuum_interpolator.initialize(mesh, poisson_solver, 0, TYPES.VACUUM);
        bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient, TYPES.BULK);
        emission.initialize(mesh);

        make_intermediate_step = true;
    } else if (make_intermediate_step) {
        surface_temperatures.calc_full_interpolation();
        make_intermediate_step = false;
    } else
        surface_temperatures.calc_interpolation();
    end_msg(t0);
}

void ProjectRunaway::make_pic_step(int& n_lost, int& n_cg, int& n_injected,
        bool full_run, bool write_files)
{
    // advance super particles
    n_lost = pic_solver.update_positions();

    // assemble and solve Poisson equation
    poisson_solver.assemble(full_run, write_files);
    n_cg = poisson_solver.solve();

    // must be after solving Poisson and before updating velocities
    vacuum_interpolator.extract_solution(poisson_solver, conf.run.field_smoother);

    // update velocities of super particles and collide them
    pic_solver.update_velocities();
    pic_solver.collide_particles();

    // update field on the surface
    surface_fields.calc_interpolation();

    // calculate field emission and inject electrons
    emission.calc_emission(conf.emission);
    n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);

    if (write_files) write_pic_results();
}

void ProjectRunaway::write_pic_results() {
    bool mode_save = MODES.WRITEFILE;
    MODES.WRITEFILE = true;

    emission.write("out/emission.dat");
    emission.write("out/emission.movie");
    pic_solver.write("out/electrons.movie");
    surface_fields.write("out/surface_fields.movie");
    surface_temperatures.write("out/surface_temperatures.movie");

//    Interpolator charge_density("elfield", "elfield_norm", "charge_density");
//    charge_density.initialize(mesh, poisson_solver, 0, TYPES.VACUUM);
//    charge_density.extract_charge_density(poisson_solver);
//    charge_density.nodes.write("out/result_E_charge.movie");

    MODES.WRITEFILE = mode_save;
    last_write_time = GLOBALS.TIME;
}

int ProjectRunaway::solve_heat(double T_ambient, double delta_time, bool full_run, int& ccg, int& hcg) {
    // Calculate field emission in case not ready from pic
    if (conf.field.mode == "laplace")
        calc_heat_emission(full_run);

    start_msg(t0, "Calculating current density");
    ch_solver.current.assemble();
    ccg = ch_solver.current.solve();
    end_msg(t0);
    check_return(ccg < 0, "Current solver did not complete normally,"
            " #CG=" + d2s(abs(ccg)) + "/" + d2s(conf.heating.n_cg));
    write_verbose_msg("#CG steps: " + d2s(ccg));

    start_msg(t0, "Calculating temperature distribution");
    ch_solver.heat.assemble(delta_time * 1.e-15); // caution!! ch_solver internal time in sec
    hcg = ch_solver.heat.solve();
    end_msg(t0);
    check_return(hcg < 0, "Heat solver did not complete normally,"
            " #CG=" + d2s(abs(hcg)) + "/" + d2s(conf.heating.n_cg));
    write_verbose_msg("#CG steps: " + d2s(hcg));

    start_msg(t0, "Extracting J & T");
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator.nodes.write("out/result_J_T.movie");
    bulk_interpolator.lintet.write("out/result_J_T.vtk");

    last_heat_time = GLOBALS.TIME;
    // TODO implement reasonable temperature limit check
    return 0;
}

void ProjectRunaway::calc_heat_emission(bool full_run) {
    static bool make_intermediate_step = false;
    start_msg(t0, "Calculating surface field & temperature");
    if (full_run) {
        surface_fields.interpolate(ch_solver);
        surface_temperatures.interpolate(ch_solver);
        bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient, TYPES.BULK); // must be after interpolation as it clears previous data
        make_intermediate_step = true;
    } else if (make_intermediate_step) {
        surface_temperatures.calc_full_interpolation();
        make_intermediate_step = false;
    } else
        surface_temperatures.calc_interpolation();
    end_msg(t0);

    surface_fields.write("out/surface_fields.movie");
    surface_temperatures.write("out/surface_temperatures.movie");

    start_msg(t0, "Calculating electron emission");
    emission.initialize(mesh, full_run);
    emission.calc_emission(conf.emission, conf.field.V0);
    end_msg(t0);

    emission.write("out/emission.movie");
}

int ProjectRunaway::write_results(bool force_write){

    if (!write_time() && !force_write) return 1;

    vacuum_interpolator.extract_solution(poisson_solver, conf.run.field_smoother);
    vacuum_interpolator.nodes.write("out/result_E_phi.movie");
    vacuum_interpolator.linhex.write("out/result_E_phi.vtk");

    if (conf.field.mode != "laplace"){
        emission.write("out/emission.movie");
        pic_solver.write("out/electrons.movie");
        surface_fields.write("out/surface_fields.movie");
        vacuum_interpolator.extract_charge_density(poisson_solver);
        vacuum_interpolator.nodes.write("out/result_E_charge.movie");
    }

    if (emission.atoms.size() > 0)
        emission.write("out/surface_emission.movie");

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

    reader.write("out/reader_err.ckx");

    if (mesh && mesh->nodes.size() > 0) {
        mesh->hexs.write("out/hexmesh_err.vtk");
        mesh->quads.write("out/quadmesh_err.vtk");
    }

    if (vacuum_interpolator.nodes.size() > 0) {
        vacuum_interpolator.nodes.write("out/result_E_phi_err.xyz");
        vacuum_interpolator.lintet.write("out/result_E_phi_err.vtk");
    }

    if (conf.heating.mode != "none" && bulk_interpolator.nodes.size() > 0) {
        bulk_interpolator.nodes.write("out/result_J_T_err.xyz");
        bulk_interpolator.lintet.write("out/result_J_T_err.vtk");
    }

    return 0;
}

int ProjectRunaway::export_data(double* data, const int n_points, const string &data_type) {
    if (n_points <= 0) return 0;

    if (fields.contains(data_type))
        return fields.export_results(n_points, data_type, data);

    if (temperatures.contains(data_type))
        return temperatures.export_results(n_points, data_type, data);

    if (forces.contains(data_type))
        return forces.export_results(n_points, data_type, data);

    if (data_type == LABELS.pair_potential_sum || data_type == LABELS.parcas_force || data_type == LABELS.charge_force)
        return forces.export_parcas(n_points, data_type, reader.get_si2parcas_box(), data);

    if (data_type == LABELS.parcas_velocity)
        return temperatures.scale_berendsen_long(data, n_points, reader.get_parcas2si_box());

    if (data_type == LABELS.atom_type) {
        require(n_points <= reader.size(), "Invalid data query size: " + d2s(n_points));
        for (int i = 0; i < n_points; ++i)
            data[i] = reader.get_marker(i);
        return reader.get_n_detached();
    }

    require(false, "Unimplemented type of export data: " + data_type);
    return 1;
}

int ProjectRunaway::interpolate(double* data, int* flag,
        const int n_points, const string &data_type, const bool near_surface,
        const double* x, const double* y, const double* z)
{
    // location of interpolation; 2-on surface, 3-in space
    int dim = 3;
    if (near_surface) dim = 2;

    if (fields.contains(data_type)) {
        fields.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        int retval = fields.interpolate_results(n_points, data_type, x, y, z, data);
        for (int i = 0; i < n_points; ++i)
            flag[i] = fields.get_marker(i) >= 0;
        return retval;
    }

    if (temperatures.contains(data_type)) {
        temperatures.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        int retval = temperatures.interpolate_results(n_points, data_type, x, y, z, data);
        for (int i = 0; i < n_points; ++i)
            flag[i] = fields.get_marker(i) >= 0;
        return retval;
    }

    require(false, "Unimplemented type of interpolation data: " + data_type);
    return 1;
}

} /* namespace femocs */
