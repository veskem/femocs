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
        fail(false), t0(0), mesh_changed(false), first_run(true),
		last_heat_time(-conf.behaviour.timestep_fs),
		last_pic_time(-conf.behaviour.timestep_fs),
		last_restart_ts(0),

        vacuum_interpolator("Elfield", "ElfieldNorm", "Potential"),
        bulk_interpolator("Rho", "Potential", "Temperature"),

        fields(&vacuum_interpolator),
        temperatures(&bulk_interpolator),
        forces(&vacuum_interpolator),

        surface_fields(&vacuum_interpolator),
        surface_temperatures(&bulk_interpolator),
        heat_transfer(&bulk_interpolator),

        phys_quantities(config.heating),
        poisson_solver(NULL, &config.field, &vacuum_interpolator.linhex),

        emission(&surface_fields, &surface_temperatures, &vacuum_interpolator),
        ch_solver(&phys_quantities, &config.heating, &emission),
        pic_solver(&poisson_solver, &emission, &vacuum_interpolator, &conf.pic, conf.behaviour.rnd_seed)
{
    poisson_solver.set_particles(pic_solver.get_particles());

    surface_fields.set_preferences(false, 2, 3);
    surface_temperatures.set_preferences(false, 2, 3);
    heat_transfer.set_preferences(false, 3, 1);

    start_msg(t0, "Reading physical quantities");
    phys_quantities.initialize_with_hc_data();
    end_msg(t0);
}

int ProjectRunaway::reinit() {
    GLOBALS.TIMESTEP++;
    conf.read_all();

    double rmsd = reader.get_rmsd();
    string rmsd_string = "inf";
    if (rmsd < 1e100) rmsd_string = d2s(rmsd);

    write_silent_msg("Running at timestep=" + d2s(GLOBALS.TIMESTEP)
            + ", time=" + d2s(GLOBALS.TIME, 2) + " fs, rmsd=" + rmsd_string);

    return rmsd < conf.geometry.distance_tol;
}

int ProjectRunaway::finalize(double tstart) {
    if (mesh_changed) reader.save_current_run_points();
    GLOBALS.TIME = GLOBALS.TIMESTEP * conf.behaviour.timestep_fs;

    write_restart("in/femocs.restart");
    write_silent_msg("Total execution time " + d2s(omp_get_wtime()-tstart, 3));

    mesh_changed = false;
    first_run = false;
    return 0;
}

int ProjectRunaway::process_failed(const string &msg) {
    write_verbose_msg(msg);
    GLOBALS.TIME = GLOBALS.TIMESTEP * conf.behaviour.timestep_fs;
    require(!first_run, "First full run failed!");

    if (conf.behaviour.n_write_log > 0)
        MODES.WRITELOG = (GLOBALS.TIMESTEP+1)%conf.behaviour.n_write_log == 0;
    return 1;
}

int ProjectRunaway::run(const int timestep, const double time) {
    double tstart = omp_get_wtime();

    //***** Build or import mesh *****

    if (reinit()) {
        write_verbose_msg("Atoms haven't moved significantly. Previous mesh will be used.");
        dense_surf.update_positions(reader);
    }

    else if (generate_mesh())
        return process_failed("Mesh generation failed!");

    if (prepare_solvers())
        return process_failed("Preparation of FEM solvers failed!");

    //***** Run FEM solvers *****

    if (run_field_solver())
        return process_failed("Running field solver in a " + conf.field.mode + " mode failed!");

    if (run_heat_solver())
        return process_failed("Running heat solver in a " + conf.heating.mode + " mode failed!");

    // it is crucial that the update of new_mesh pointer AND mesh of bulk_interpolator
    // happens within time interval where no error can occur
    update_mesh_pointers();

    //***** Prepare for data export and next run *****

    if (prepare_export())
        return process_failed("Interpolating solution on atoms failed!");

    return finalize(tstart);
}

int ProjectRunaway::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum) {
    start_msg(t0, "Extracting surface");
    reader.extract(dense_surf, TYPES.SURFACE);
    end_msg(t0);
    dense_surf.write("out/surface_dense.xyz");

    if (first_run) {
        start_msg(t0, "Extending surface");
        dense_surf.extend(extended_surf, conf);
        end_msg(t0);
        extended_surf.write("out/surface_extension.xyz");
    }

    start_msg(t0, "Coarsening surface");
    dense_surf.generate_boundary_nodes(bulk, coarse_surf, vacuum, extended_surf, conf, first_run);
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

    new_mesh->nodes.write("out/hexmesh_nodes.vtk");
    new_mesh->tris.write("out/trimesh.vtk");
    new_mesh->quads.write("out/quadmesh.vtk");
    new_mesh->tets.write("out/tetmesh.vtk");
    new_mesh->hexs.write("out/hexmesh.vtk");
    new_mesh->write_separate("out/hexmesh_bulk.vtks", TYPES.BULK);

    // update mesh pointers
    mesh = new_mesh;
    mesh_changed = true;

    write_verbose_msg(mesh->to_str());
    return 0;
}

int ProjectRunaway::import_mesh() {
    start_msg(t0, "Importing vacuum mesh to Deal.II");
    fail = !poisson_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);

    if (conf.field.mode != "laplace" || conf.heating.mode != "none") {
        start_msg(t0, "Importing bulk mesh to Deal.II");
        fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexs.export_bulk());
        check_return(fail, "Importing bulk mesh to Deal.II failed!");
        end_msg(t0);
    }

    return 0;
}

void ProjectRunaway::update_mesh_pointers() {
    static bool odd_run = true;
    if (mesh_changed) {
        if (odd_run) new_mesh = &mesh2;
        else new_mesh = &mesh1;
        odd_run = !odd_run;
    }
}

void ProjectRunaway::calc_surf_temperatures() {
    static bool make_intermediate_step = false;
    if (mesh_changed) {
        surface_temperatures.interpolate(ch_solver);
        make_intermediate_step = true;
    } else if (make_intermediate_step) {
        surface_temperatures.calc_full_interpolation();
        make_intermediate_step = false;
    } else
        surface_temperatures.calc_interpolation();
}

int ProjectRunaway::prepare_solvers() {
    // import Tetgen mesh to Deal.II
    if (mesh_changed) {
        int error = import_mesh();
        if (error) return 1;
    }

    // halt in case of NO field emission
    if (conf.field.mode == "laplace" && conf.heating.mode == "none")
        return 0;

    // initialize the calculation of field emission
    if (mesh_changed) {
        // setup ch_solver here as it must be done before transferring previous heat values into new mesh,
        // which in turn must be done before re-initializing bulk_interpolator
        start_msg(t0, "Setup current & heat solvers");
        ch_solver.setup(conf.heating.t_ambient);
        end_msg(t0);

        ch_solver.export_surface_centroids(surface_fields);
        emission.initialize(mesh);

        // in case of first run, give initial values to the temperature interpolator,
        // otherwise transfer existing temperatures to new mesh
        if (bulk_interpolator.nodes.size() == 0) {
            bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient, TYPES.BULK);
        } else {
            start_msg(t0, "Transferring old temperatures to new mesh");
            heat_transfer.interpolate_dofs(ch_solver);
            end_msg(t0);
        }
    }

    start_msg(t0, "Calculating surface temperatures");
    calc_surf_temperatures();
    end_msg(t0);

    return 0;
}

int ProjectRunaway::prepare_export() {
    if (dense_surf.size() == 0) {
        write_silent_msg("No atoms for solution interpolation!");
        return 0;
    }

    start_msg(t0, "Interpolating E & phi");
    if (mesh_changed) {
        fields.set_preferences(false, 2, conf.behaviour.interpolation_rank);
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
            temperatures.precalc_berendsen();
        } else {
            temperatures.update_positions(reader);
        }
        end_msg(t0);
        // writing temperatures will be done during Berendsen scaling
        // TODO implement reasonable temperature limit check
    }

    int retval = 0;
    if (conf.force.mode != "none") {
        if (mesh_changed)
            retval = solve_force();
        else {
            start_msg(t0, "Recalculating forces");
            forces.update_positions(dense_surf);
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
        return solve_pic(GLOBALS.TIME - last_pic_time, mesh_changed);

    else if (mesh_changed)
        return solve_laplace(conf.field.E0, conf.field.V0);

    return 0;
}

int ProjectRunaway::run_heat_solver() {
    int ccg, hcg;

    double delta_time = GLOBALS.TIME - last_heat_time;
    bool b1 = delta_time >= conf.heating.delta_time;
    if (conf.heating.mode == "transient" && (mesh_changed || b1))
        return solve_heat(conf.heating.t_ambient, delta_time, mesh_changed, ccg, hcg);

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
    poisson_solver.assemble(true);
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
    pic_solver.set_params(dt_pic, mesh->nodes.stat);
    GLOBALS.TIME = GLOBALS.TIMESTEP * conf.behaviour.timestep_fs - advance_time;
    last_pic_time = GLOBALS.TIME - conf.behaviour.timestep_fs;

    if (full_run) {
        start_msg(t0, "Initializing Poisson solver");
        poisson_solver.setup(-conf.field.E0, conf.field.V0);
        vacuum_interpolator.initialize(mesh, poisson_solver, 0, TYPES.VACUUM);
        end_msg(t0);
        write_verbose_msg(poisson_solver.to_str());
    }

    int n_lost, n_cg, n_injected, error;
    start_msg(t0, "=== Running PIC for delta time = "
            + d2s(n_pic_steps) + "*" + d2s(dt_pic)+" = " + d2s(advance_time) + " fs\n");

    for (int i = 0; i < n_pic_steps; ++i) {
        error = make_pic_step(n_lost, n_cg, n_injected, full_run && i==0);
        if (error) return 1;

        if (MODES.VERBOSE) {
            printf("  #CG=%d, Fmax=%.3f V/A, Itot=%.3e A, M/A=%.3f, #el_inj|del|tot=%d|%d|%d\n",
                    n_cg, emission.global_data.Fmax, emission.global_data.I_tot, fields.get_beta(),
                    n_injected, n_lost, pic_solver.get_n_electrons());
        }

        GLOBALS.TIME += dt_pic;
        last_pic_time += dt_pic;
    }
    end_msg(t0);

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintet.write("out/result_E_phi.vtk");

    return 0;
}

int ProjectRunaway::make_pic_step(int& n_lost, int& n_cg, int& n_injected, bool full_run) {
    // advance super particles
    n_lost = pic_solver.update_positions();

    // assemble and solve Poisson equation
    poisson_solver.assemble(full_run);
    n_cg = poisson_solver.solve();

    // must be after solving Poisson and before updating velocities
    vacuum_interpolator.extract_solution(poisson_solver, conf.run.field_smoother);
    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions(), false),
            "Field enhancement is out of limits!");

    // update velocities of super particles and collide them
    pic_solver.update_velocities();
    pic_solver.collide_particles();

    // update field on the surface
    surface_fields.calc_interpolation();

    // calculate field emission and inject electrons
    emission.calc_emission(conf.emission, conf.field.V0);
    n_injected = abs(pic_solver.inject_electrons(conf.pic.fractional_push));
    check_return(n_injected > conf.pic.max_injected,
            "Too many injected SP-s, " + d2s(n_injected) + ". Check the SP weight!")

    emission.write("out/emission.dat", FileIO::no_update);
    emission.write("out/emission.movie");
    pic_solver.write("out/electrons.movie");
    surface_fields.write("out/surface_fields.movie");
    surface_temperatures.write("out/surface_temperatures.movie");

    // NB! If charge_density really needed, move it into header
//    Interpolator charge_density("Elfield", "ElfieldNorm", "ChargeDensity");
//    if (charge_density.nodes.write_time()) {
//        charge_density.initialize(mesh, poisson_solver, 0, TYPES.VACUUM);
//        charge_density.extract_charge_density(poisson_solver);
//        charge_density.nodes.write("out/result_E_charge.movie");
//    }

    return 0;
}

int ProjectRunaway::solve_heat(double T_ambient, double delta_time, bool full_run, int& ccg, int& hcg) {
    // Calculate field emission in case not ready from PIC
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

    ch_solver.write("out/ch_solver.movie");
    write_verbose_msg("#CG steps: " + d2s(hcg));

    start_msg(t0, "Extracting J & T");
    bulk_interpolator.initialize(mesh, ch_solver, conf.heating.t_ambient, TYPES.BULK);
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    bulk_interpolator.nodes.write("out/result_J_T.xyz");
    bulk_interpolator.lintet.write("out/result_J_T.vtk");

    last_heat_time = GLOBALS.TIME;
    // TODO implement reasonable temperature limit check
    return 0;
}

void ProjectRunaway::calc_heat_emission(bool full_run) {
    if (full_run) {
        start_msg(t0, "Calculating surface fields");
        surface_fields.calc_interpolation();
        end_msg(t0);
        surface_fields.write("out/surface_fields.movie");
    }

    start_msg(t0, "Calculating electron emission");
    emission.calc_emission(conf.emission, conf.field.V0);
    end_msg(t0);
    emission.write("out/emission.movie");
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
        return temperatures.scale_berendsen(data, n_points, reader.get_parcas2si_box(), conf);

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

int ProjectRunaway::restart(const string &path_to_file) {
    start_msg(t0, "Reading restart file from " + path_to_file);
    fail = new_mesh->read(path_to_file, "");
    bulk_interpolator.initialize(new_mesh, TYPES.BULK);
    bulk_interpolator.nodes.read(path_to_file, 1);
    pic_solver.read(path_to_file);

    last_pic_time = GLOBALS.TIME - conf.behaviour.timestep_fs;
    last_heat_time = GLOBALS.TIME - conf.behaviour.timestep_fs;
    last_restart_ts = GLOBALS.TIMESTEP;
    end_msg(t0);

    write_verbose_msg(new_mesh->to_str());

    require(new_mesh->nodes.size() == bulk_interpolator.nodes.size(),
            "Mismatch between mesh and interpolator sizes: " +
            d2s(new_mesh->nodes.size()) + " vs " + d2s(bulk_interpolator.nodes.size()));

    mesh = new_mesh;
    mesh_changed = true;
    update_mesh_pointers();

    return 0;
}

void ProjectRunaway::write_restart(const string &path_to_file) {
    if (conf.behaviour.n_write_restart <= 0 ||
             (GLOBALS.TIMESTEP - last_restart_ts) < conf.behaviour.n_write_restart)
        return;

    unsigned int flags = FileIO::force | FileIO::no_update | FileIO::append;
    mesh->write(path_to_file, FileIO::force | FileIO::no_update);
    bulk_interpolator.nodes.write(path_to_file, flags);
    pic_solver.write(path_to_file, flags);

    last_restart_ts = GLOBALS.TIMESTEP;
}

} /* namespace femocs */
