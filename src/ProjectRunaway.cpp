/*
 * ProjectRunaway.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include <omp.h>
#include <algorithm>
#include <sstream>
#include <cmath>

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
        emission(&surface_fields, &surface_temperatures, &vacuum_interpolator),

        phys_quantities(config.heating),
        poisson_solver(NULL, &config.field),
        ch_solver(&phys_quantities, &config.heating),
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

double ProjectRunaway::max_field() {
    double max_field = 0;
    for (Solution s : *vacuum_interpolator.nodes.get_solutions())
        max_field = max(max_field, s.norm);
    return max_field;
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

    //************************* MESHING ************************
    bool skip_meshing = reinit(tstep);
    if (skip_meshing)
        write_verbose_msg(stream.str());
    else {
        if (generate_mesh()) {
            force_output();
            check_return(true, "Mesh generation failed!");
        }
        if (prepare_solvers()) {
            force_output();
            check_return(true, "Preparation of FEM solvers failed!");
        }
    }

    //****************** RUNNING FEM SOLVERS *******************
    if (conf.field.solver == "poisson") {
        if (solve_pic(elfield)) {
            force_output();
            check_return(true, "Solving PIC failed!");
        }
    }

    if (!skip_meshing && conf.field.solver == "laplace") {
        if (solve_laplace(elfield)) {
            force_output();
            check_return(true, "Solving Laplace equation failed!");
        }
    }

    if (!skip_meshing && solve_heat(conf.heating.t_ambient)) {
        force_output();
        check_return(true, "Solving heat & continuity equation failed!");
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

        surface_fields.set_preferences(false, 2, conf.behaviour.interpolation_rank);
        surface_temperatures.set_preferences(false, 2, conf.behaviour.interpolation_rank);
    }

    return 0;
}

int ProjectRunaway::solve_laplace(const double E0) {
    bool first_time = true;
    conf.field.E0 = E0;       // reset long-range electric field

    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(E0, conf.tolerance.field_min, conf.tolerance.field_max, conf.geometry.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Initializing Poisson solver...");
    poisson_solver.setup(-E0);
    poisson_solver.assemble_laplace(first_time);
    end_msg(t0);

    stringstream ss; ss << poisson_solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Poisson solver...");
    int ncg = poisson_solver.solve();
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    vacuum_interpolator.extract_solution(poisson_solver);
    end_msg(t0);

    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()), "Field enhancement is out of limits!");

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi_linear.vtk");
    vacuum_interpolator.quadtets.write("out/result_E_phi_quad.vtk");

    return fail;
}

int ProjectRunaway::solve_pic(const double E0) {
    const double dt_main = max(delta_t_MD * 1.e15, conf.pic.total_time);

    int time_subcycle = ceil(dt_main / conf.pic.dt_max); // dt_main = delta_t_MD converted to [fs]
    double dt_pic = dt_main/time_subcycle;

    start_msg(t0, "=== Interpolating elfield on faces...");
    surface_fields.interpolate(ch_solver);
    end_msg(t0);
    
    surface_fields.write("out/surface_field.xyz");
    
    // Temperature do not need to be recalculated, as thermal timestep is >> PIC timestep
    start_msg(t0, "=== Interpolating J & T on faces...");
    surface_temperatures.interpolate(ch_solver);
    end_msg(t0);
    
    surface_temperatures.write("out/surface_temperature.xyz");

    start_msg(t0, "=== Initializing Poisson solver...");
    poisson_solver.setup(-E0, conf.field.V0);
    pic_solver.set_params(conf.field, conf.pic, dt_pic, mesh->nodes.stat);
    emission.initialize(mesh);
    end_msg(t0);
    
    start_msg(t0, "=== Running PIC...\n");
    int n_lost, n_injected, n_cg_steps;
    
    for (int i = 0; i < time_subcycle; i++) {
        n_lost = pic_solver.update_positions();
        n_cg_steps = pic_solver.run_cycle(i == 0);

        vacuum_interpolator.extract_solution(poisson_solver);
        surface_fields.calc_interpolation();
        emission.calc_emission(conf.emission, conf.field.V0);

        n_injected = pic_solver.inject_electrons(conf.pic.fractional_push);
        
        if (MODES.VERBOSE)
            printf("  #CG steps=%d, max field=%.3f, #injected|deleted electrons=%d|%d\n",
                    n_cg_steps, max_field(), n_injected, n_lost);
    }
    
    end_msg(t0);

    pic_solver.write("out/electrons.movie", 0);

    return 0;

    //7. Save ions and neutrals that are inbound on the MD domain somewhere where the MD can find them
    // TODO LATER
    //8. Give the heat- and current fluxes to the temperature solver.
    // TODO LATER
}

int ProjectRunaway::solve_heat(const double T_ambient) {
    if(conf.heating.mode == "transient") {
        double delta_time = delta_t_MD * (timestep - last_full_timestep); //in sec
        return solve_transient_heat(T_ambient, delta_time);
    }
    else if (conf.heating.mode == "converge") {
        return solve_converge_heat(T_ambient);
    }

    return 0;
}

int ProjectRunaway::solve_transient_heat(const double T_ambient, const double delta_time) {
    double multiplier = 1.;

    // Interpolate elfield on face centroids
    surface_fields.interpolate(ch_solver);

    // Interpolate J & T on face centroids
    surface_temperatures.interpolate(ch_solver);

    // Calculate field emission
    emission.initialize(mesh);
    emission.calc_emission(conf.emission, conf.field.V0);
    emission.export_emission(ch_solver);

    start_msg(t0, "=== Setup transient J & T solver...");
    ch_solver.setup(T_ambient);
    end_msg(t0);

    start_msg(t0, "=== Calculating current density...");
    ch_solver.current.assemble();
    unsigned int ccg = ch_solver.current.solve();
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(ccg));

    start_msg(t0, "=== Calculating temperature distribution...");
    ch_solver.heat.assemble(delta_time);
    unsigned int hcg = ch_solver.heat.solve();
    end_msg(t0);
    write_verbose_msg("# CG steps: " + to_string(hcg));

    start_msg(t0, "=== Extracting J & T...");
    bulk_interpolator.extract_solution(ch_solver);
    end_msg(t0);

    write_verbose_msg("Current and heat advanced for " + to_string(delta_time));

    bulk_interpolator.nodes.write("out/result_J_T.xyz");
    bulk_interpolator.lintets.write("out/result_J_T.vtk");

    return 0;
}

int ProjectRunaway::solve_converge_heat(const double T_ambient) {
    start_msg(t0, "=== Interpolating elfield on faces...");
    surface_fields.interpolate(ch_solver);
    end_msg(t0);

    surface_fields.write("out/surface_field.xyz");

    start_msg(t0, "=== Setup transient J & T solver...");
    ch_solver.setup(T_ambient);
    end_msg(t0);

    emission.initialize(mesh);

    double current_time = 0.;
    double delta_time = 1.e-12; //in seconds!!
    double multiplier = 1.;

    start_msg(t0, "=== Running converge J & T solver...\n");
    int converge_steps = 0;
    for (; converge_steps < 1000; ++converge_steps) {
        double tstart = omp_get_wtime();

        // Interpolate J & T on face centroids
        surface_temperatures.interpolate(ch_solver);

        // Calculating field emission
        multiplier = emission.calc_emission(multiplier, conf.emission, conf.field.V0);
        emission.export_emission(ch_solver);

        // Calculate current density
        ch_solver.current.assemble();
        unsigned int ccg = ch_solver.current.solve(); // ccg == nr of CG iterations

        // Calculate temperature distribution
        ch_solver.heat.assemble(delta_time);
        unsigned int hcg = ch_solver.heat.solve(); // hcg == nr of CG iterations

        // Extract J & T
        bulk_interpolator.extract_solution(ch_solver);

        if (MODES.VERBOSE) {
            double max_T = ch_solver.heat.max_solution();
            double time = omp_get_wtime() - tstart;
            printf("  i=%d, t=%.2fps, dt=%.2fps, ccg=%d, hcg=%d, time=%.3fs, Tmax=%.2fK\n",
                    converge_steps, current_time * 1.e12, delta_time * 1.e12, ccg, hcg, time, max_T);
        }
        current_time += delta_time;

        if (max(hcg, ccg) < 120 || hcg < 30)
            delta_time *= 1.25;
        else if (max(hcg, ccg) > 150)
            delta_time /= 1.25;

        if (max(hcg, ccg) < 10) break;
    }
    end_msg(t0);
    if (converge_steps >= 1000)
        write_silent_msg("WARNING: Heat equation did not converge after 1000 steps!");

    surface_temperatures.write("out/surface_temperature.xyz");
    emission.write("out/surface_emission.xyz");
    bulk_interpolator.nodes.write("out/result_J_T.xyz");

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
    if (cmd == LABELS.elfield.second || cmd == LABELS.elfield_norm.second || cmd == LABELS.potential.second)
        return fields.export_results(n_points, cmd, false, data);

    require(false, "Unimplemented type of export data: " + cmd);
    return 0;
}

int ProjectRunaway::interpolate_results(const int n_points, const string &cmd, const bool surface,
        const double* x, const double* y, const double* z, double* data, int* flag) {

    if (cmd == LABELS.elfield.second || cmd == LABELS.elfield_norm.second || cmd == LABELS.potential.second) {
        fields.set_preferences(true, 3, conf.behaviour.interpolation_rank);
        return fields.interpolate_results(n_points, cmd, x, y, z, data);
    }

    require(false, "Unimplemented type of export data: " + cmd);
    return 0;
}

} /* namespace femocs */
