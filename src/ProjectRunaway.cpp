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
        poisson_solver(NULL, &config.field),
        ch_solver(&phys_quantities, &config.heating),

        emission(&surface_fields, &surface_temperatures, &poisson_solver, &vacuum_interpolator),
        pic_solver(&poisson_solver, &ch_solver, &vacuum_interpolator, &emission)
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

int ProjectRunaway::reinit(const int tstep) {
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

    bool skip_meshing = reader.rms_distance < conf.tolerance.distance;
    prev_skip_calculations = skip_meshing;
    return skip_meshing;
}

int ProjectRunaway::finalize() {
    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.tolerance.distance);
    end_msg(t0);
    last_full_timestep = timestep;
    return 0;
}

int ProjectRunaway::run(const int timestep) {
    return run(conf.field.E0, timestep);
}

int ProjectRunaway::run(const double elfield, const int tstep) {
    stringstream stream;
    stream << fixed << setprecision(3);

    double tstart = omp_get_wtime();

    stream.str("");
    stream << "Atoms haven't moved significantly, " << reader.rms_distance
            << " < " << conf.tolerance.distance << "! Previous mesh will be used!";

    //***** Build mesh *****
    bool skip_meshing = reinit(tstep);
    if (skip_meshing)
        write_verbose_msg(stream.str());
    else {
        new_mesh_exists = true;
        if (generate_mesh()) {
            force_output();
            write_silent_msg("Mesh generation failed!");
            new_mesh_exists = false;
        }
        if (prepare_solvers()) {
            force_output();
            write_silent_msg("Preparation of FEM solvers failed!");
            new_mesh_exists = false;
        }
    }
    check_return(time == 0 && !new_mesh_exists, "first meshing failed! FEMOCS terminating...");

    //***** Run heat convergence calculation with constant mesh
    if (conf.heating.mode == "converge") return run_heat_converge(conf.heating.t_ambient);

    //***** Run FEM solvers *****

    if (conf.field.solver == "poisson") { // solver Poisson equation - run PIC
        if (solve_pic(conf.behaviour.total_time)) {
            force_output();
            check_return(true, "Solving PIC failed!");
        }
    }

    if (new_mesh_exists && conf.field.solver == "laplace") {
        if (solve_laplace(elfield)) {
            force_output();
            check_return(true, "Solving Laplace equation failed!");
        }
    }

    int ccg, hcg;
    if (!skip_meshing && solve_transient_heat(conf.heating.t_ambient, time - last_heat_time, ccg, hcg)) {
        force_output();
        check_return(true, "Solving heat & continuity equation failed!");
    }

    // In case no transient simulations (time has stayed the same) advance the time
    if (!conf.pic.run_pic && conf.heating.mode == "none")
        time += conf.behaviour.total_time;

    //***** Prepare for data export and next run *****
    if (prepare_export()) {
        force_output();
        check_return(true, "Extracting solution on atoms failed!");
    }

    if (is_write_time()){
        start_msg(t0, "=== Writing to files...");
        write();
        end_msg(t0);
    }

    finalize();

    stream.str("");
    stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());

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
    check_return(err_code, "Triangulation failed with error code " + to_string(err_code));
    end_msg(t0);

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = new_mesh->mark_mesh();
    check_return(fail, "Mesh marking failed!");
    end_msg(t0);

    start_msg(t0, "=== Generating surface faces...");
    err_code = new_mesh->generate_surface(reader.sizes, "rQB", "rQnn");
    end_msg(t0);
    check_return(err_code, "Generation of surface faces failed with error code " + to_string(err_code));

    if (conf.smoothing.algorithm != "none" && conf.smoothing.n_steps > 0) {
        start_msg(t0, "=== Smoothing surface faces...");
        new_mesh->smoothen(conf.smoothing.n_steps, conf.smoothing.lambda_mesh, conf.smoothing.mu_mesh, conf.smoothing.algorithm);
        end_msg(t0);
    }

    // has to be separate to ensure that data is for sure calculated
    new_mesh->faces.calc_norms_and_areas();

    new_mesh->nodes.write("out/tetmesh_nodes.vtk");
    new_mesh->faces.write("out/trimesh.vtk");
    new_mesh->elems.write("out/tetmesh.vtk");

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
    new_mesh->hexahedra.write("out/hexmesh.vtk");
    new_mesh->write_separate("out/hexmesh_bulk" + timestep_string + ".vtk", TYPES.BULK);
    new_mesh->faces.write("out/hexmesh_faces.vtk");

    // update mesh pointers
    static bool odd_run = true;

    mesh = new_mesh;
    if (odd_run) new_mesh = &mesh2;
    else new_mesh = &mesh1;

    odd_run = !odd_run;

    stringstream ss; ss << *mesh;
    write_verbose_msg(ss.str());

    return 0;
}

int ProjectRunaway::prepare_solvers() {
    start_msg(t0, "=== Importing mesh to Poisson solver...");
    fail = !poisson_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexahedra.export_vacuum());
    check_return(fail, "Importing vacuum mesh to Deal.II failed!");
    end_msg(t0);

    vacuum_interpolator.initialize(mesh);
    vacuum_interpolator.lintets.narrow_search_to(TYPES.VACUUM);

    if (conf.field.solver == "poisson" || conf.heating.mode != "none") {
        start_msg(t0, "=== Importing mesh to J & T solver...");
        fail = !ch_solver.import_mesh(mesh->nodes.export_dealii(), mesh->hexahedra.export_bulk());
        check_return(fail, "Importing bulk mesh to Deal.II failed!");
        end_msg(t0);

        bulk_interpolator.initialize(mesh, conf.heating.t_ambient);
        bulk_interpolator.lintets.narrow_search_to(TYPES.BULK);

        surface_fields.set_preferences(false, 2, 3);
        surface_temperatures.set_preferences(false, 2, 3);
    }

    return 0;
}

int ProjectRunaway::prepare_export() {
    start_msg(t0, "=== Interpolating E and phi...");
    fields.set_preferences(true, 2, conf.behaviour.interpolation_rank);
    fields.interpolate(dense_surf);
    end_msg(t0);

    fields.write("out/fields.movie");
    check_return(fields.check_limits(), "Field enhancement is out of limits!");

    if (conf.heating.mode != "none") {
        start_msg(t0, "=== Interpolating J & T...");
        temperatures.set_preferences(true, 3, conf.behaviour.interpolation_rank);
        temperatures.interpolate(reader);
        end_msg(t0);

        temperatures.write("out/temperatures.movie");
        // TODO implement reasonable temperature limit check
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
        forces.distribute_charges(fields, face_charges, 0, conf.smoothing.beta_charge);
        end_msg(t0);

        start_msg(t0, "=== Generating Voronoi cells...");
        VoronoiMesh voro_mesh;
        int err_code;
        err_code = forces.calc_voronois(voro_mesh, atom2face, conf.geometry.radius, conf.geometry.latconst, "10.0");
        end_msg(t0);

        check_return(err_code, "Generation of Voronoi cells failed with error code " + to_string(err_code));
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

int ProjectRunaway::solve_laplace(const double E0) {
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(E0, conf.tolerance.field_min, conf.tolerance.field_max, conf.geometry.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Initializing Laplace solver...");
    poisson_solver.setup(-E0);
    poisson_solver.assemble_laplace(true);
    end_msg(t0);

    stringstream ss; ss << poisson_solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Laplace solver...");
    int ncg = poisson_solver.solve();
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    vacuum_interpolator.extract_solution(poisson_solver);
    end_msg(t0);

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi.vtk");
    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()), "Field enhancement is out of limits!");

    return 0;
}

int ProjectRunaway::solve_pic(double advance_time) {
    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(conf.field.E0, conf.tolerance.field_min, conf.tolerance.field_max,
            conf.geometry.radius, dense_surf.sizes.zbox);

    int time_subcycle = ceil(advance_time / conf.pic.dt_max); // dt_main = delta_t_MD converted to [fs]
    double dt_pic = advance_time / time_subcycle;


    if (new_mesh_exists){
        start_msg(t0, "=== Initializing Poisson solver...");
        poisson_solver.setup(conf.field.E0, conf.field.V0);
        end_msg(t0);
    }
    pic_solver.set_params(conf.field, conf.pic, dt_pic, mesh->nodes.stat);
    
    if (conf.pic.run_pic)
        start_msg(t0, "=== Running PIC...\n");
    else
        start_msg(t0, "=== Solving Poisson equation...\n");
    
    for (int i = 0; i < time_subcycle; i++) {
        int n_lost = pic_solver.update_positions();
        int n_cg_steps = pic_solver.run_cycle(new_mesh_exists, is_write_time());

        // If not PIC, advance time by whole timestep and return
        if (!conf.pic.run_pic)
            break;

        // Set up emission input
        if (i == 0) {
            if (new_mesh_exists) surface_fields.interpolate(ch_solver);
            surface_temperatures.interpolate(ch_solver);
            emission.initialize(mesh);
        }

        //calculate emission and inject electrons
        emission.calc_emission(conf.emission, conf.field.V0);
        int n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);
        
        if (MODES.VERBOSE)
            printf("t= %f fs, #CG =%d, Fmax=%.3f V/A, Itot=%.3e A #e|inj|del|tot=%d|%d|%d\n",
                    time, n_cg_steps, emission.global_data.Fmax, emission.global_data.I_tot, n_injected,
                    n_lost, pic_solver.get_n_electrons());

        time += dt_pic;

        if (is_write_time()){
            start_msg(t0, "=== Writing to files...");
            write();
            end_msg(t0);
        }

    }
    
    end_msg(t0);
    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()),
            "Field enhancement is out of limits!");

    return 0;

    //7. Save ions and neutrals that are inbound on the MD domain somewhere where the MD can find them
    // TODO LATER
    //8. Give the heat- and current fluxes to the temperature solver.
    // TODO LATER
}


int ProjectRunaway::solve_transient_heat(const double T_ambient, const double delta_time, int& ccg, int& hcg) {

    if (new_mesh_exists)
        ch_solver.setup(T_ambient);

    // Calculate field emission in case not ready from pic
    if(!conf.pic.run_pic){
        start_msg(t0, "=== Calculating electron emission...");
        if (new_mesh_exists) surface_fields.interpolate(ch_solver);
        surface_temperatures.interpolate(ch_solver);
        emission.initialize(mesh);
        emission.calc_emission(conf.emission, conf.field.V0);
        emission.export_emission(ch_solver);
        end_msg(t0);
    }

    start_msg(t0, "=== Calculating current density...");
    ch_solver.current.assemble();
    ccg = ch_solver.current.solve();
    end_msg(t0);
    write_verbose_msg("#CG steps: " + to_string(ccg));

    start_msg(t0, "=== Calculating temperature distribution...");
    ch_solver.heat.assemble(delta_time * 1.e-15); // caution!! ch_solver internal time in sec
    hcg = ch_solver.heat.solve();
    end_msg(t0);
    write_verbose_msg("#CG steps: " + to_string(hcg));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    if(!conf.pic.run_pic) time += delta_time;

    if (is_write_time()){
        start_msg(t0, "=== Writing to files...");
        write();
        end_msg(t0);
    }

    last_heat_time = time;

    return 0;
}

int ProjectRunaway::run_heat_converge(const double T_ambient) {

    double delta_time = 1;

    start_msg(t0, "=== Running heat convergence...\n");
    int ccg, hcg;
    int converge_steps = 0;
    for (; converge_steps < 1000; ++converge_steps) {

        //calculate field - advance pic for the given timestep
        if (conf.pic.run_pic || converge_steps == 0){
            if (solve_pic(delta_time)) {
                force_output();
                check_return(true, "Solving PIC failed!");
            }
        }

        // advance heat and current system by delta_time
        if (solve_transient_heat(conf.heating.t_ambient, delta_time, ccg, hcg)) {
            force_output();
            check_return(true, "Solving heat & continuity equation failed!");
        }

        new_mesh_exists = false; //running convergence with the same mesh always

        if (MODES.VERBOSE) {
            double max_T = ch_solver.heat.max_solution();
            printf("t=%.2e ps, dt=%.2f ps, Tmax=%.2fK\n",
                    time * 1.e-3, delta_time * 1.e-3, ccg, hcg, max_T);
        }

        if (!conf.pic.run_pic)
            time += delta_time;

        if (hcg < (ccg - 10)) // if heat changed too little
            delta_time *= 1.25;
        else if (hcg > (ccg + 10)) // if heat changed too much
            delta_time /= 1.25;

        if (is_write_time()){
            start_msg(t0, "=== Writing to files...");
            write();
            end_msg(t0);
        }


        if (max(hcg, ccg) < 10) break;
    }
    end_msg(t0);

    if (converge_steps >= 1000)
        write_silent_msg("WARNING: Heat equation did not converge after 1000 steps!");


    return 0;
}

int ProjectRunaway::force_output() {
    if (conf.behaviour.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    mesh->hexahedra.write("out/hexmesh_err.vtk");
    mesh->elems.write("out/tetmesh_err.vtk");
    mesh->faces.write("out/trimesh_err.vtk");

    vacuum_interpolator.nodes.write("out/result_E_phi_err.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi_err.vtk");

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
        return fields.interpolate_results(n_points, cmd, x, y, z, data);
    }

    if (cmd == LABELS.force.second || cmd == LABELS.force_norm.second || cmd == LABELS.charge.second) {
        forces.set_preferences(false, dim, conf.behaviour.interpolation_rank);
        return forces.interpolate_results(n_points, cmd, x, y, z, data);
    }

    require(false, "Unimplemented type of interpolation data: " + cmd);
    return 1;
}


int ProjectRunaway::write(){
    // write the field
    vacuum_interpolator.extract_solution(poisson_solver);
    vacuum_interpolator.nodes.write("out/result_E_phi.movie");
    vacuum_interpolator.lintets.write("out/result_E_phi.vtk");

    // write PIC (particles and charge density)
    if (conf.pic.run_pic){
        pic_solver.write("out/electrons.movie", 0);
        vacuum_interpolator.extract_charge_density(poisson_solver);
        vacuum_interpolator.nodes.write("out/result_E_charge.movie");
    }

    //write emission
    if (emission.atoms.size() > 0){
        emission.export_emission(ch_solver);
        emission.write_data("out/emission_data.dat", time);
        emission.write("out/surface_emission.movie");
    }

    //write heat
    if (conf.heating.mode != "none"){
        bulk_interpolator.nodes.write("out/result_J_T.movie");
        bulk_interpolator.lintets.write("out/result_J_T.vtk");
    }

    last_write_time = time;

    return 0;
}

} /* namespace femocs */
