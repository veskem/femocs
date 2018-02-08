/*
 * ProjectNanotip.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include "GeneralProject.h"
#include "ProjectNanotip.h"
#include "Coarseners.h"
#include "Macros.h"
#include "Tethex.h"
#include "VoronoiMesh.h"
#include "TetgenMesh.h"

#include <omp.h>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;
namespace femocs {

ProjectNanotip::ProjectNanotip(AtomReader &a, Config &c) :
        GeneralProject(a, c),
        skip_meshing(false), fail(false), t0(0), timestep(-1), last_full_timestep(0)
{
    fields.set_interpolator(&vacuum_interpolator);
}

// Write all the available data to file for debugging purposes
int ProjectNanotip::force_output() {
    if (conf.behaviour.n_writefile <= 0) return 1;

    MODES.WRITEFILE = true;

    reader.write("out/reader.xyz");
    mesh->hexahedra.write("out/hexmesh_err.vtk");
    mesh->elems.write("out/tetmesh_err.vtk");
    mesh->faces.write("out/trimesh_err.vtk");

    vacuum_interpolator.nodes.write("out/result_E_phi_err.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi_vacuum_err.vtk");

    return 0;
}

int ProjectNanotip::run(const int timestep) {
    return run(conf.field.E0, timestep);
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int ProjectNanotip::run(const double elfield, const int tstep) {
    stringstream stream;
    stream << fixed << setprecision(3);

    double tstart = omp_get_wtime();

    stream.str("");
    stream << "Atoms haven't moved significantly, " << reader.rms_distance
            << " < " << conf.tolerance.distance << "! Previous mesh will be used!";

    //******************** MESHING *****************************
    if (reinit(tstep)) { // reinit and check skip_meshing
        write_verbose_msg(stream.str());
    } else {
        generate_mesh();
    }

    //****************** RUNNING Field - PIC calculation ********
    skip_meshing = true;

    if (solve_laplace(elfield)) {
        force_output();
        check_return(true, "Solving Laplace equation failed!");
    }

    finalize();

    stream.str("");
    stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());

    return 0;
}

// Determine whether atoms have moved significantly and whether to enable file writing
int ProjectNanotip::reinit(const int tstep) {
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

    prev_skip_calculations = skip_meshing;
    return skip_meshing;
}

// Store the imported atom coordinates and set flag that enables exporters
int ProjectNanotip::finalize() {
    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.tolerance.distance);
    end_msg(t0);
    skip_meshing = false;
    last_full_timestep = timestep;
    return 0;
}

// Generate boundary nodes for mesh
int ProjectNanotip::generate_boundary_nodes(Surface& bulk, Surface& coarse_surf, Surface& vacuum) {
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

// Generate bulk and vacuum meshes
int ProjectNanotip::generate_mesh() {
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

// Solve Laplace equation
int ProjectNanotip::solve_laplace(const double E0) {
    conf.field.E0 = E0;       // reset long-range electric field

    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    fields.set_check_params(E0, conf.tolerance.field_min, conf.tolerance.field_max, conf.geometry.radius, dense_surf.sizes.zbox);

    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !laplace_solver.import_mesh_directly(mesh->nodes.export_dealii(),
            mesh->hexahedra.export_vacuum());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Initializing Laplace solver...");
    laplace_solver.setup_system(true);
    laplace_solver.assemble_system(-E0);
    end_msg(t0);

    stringstream ss; ss << laplace_solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Laplace solver...");
    int ncg = laplace_solver.solve(conf.field.n_phi, conf.field.phi_error, true, conf.field.ssor_param);
    cout << "CG iterations = " << ncg;
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    vacuum_interpolator.initialize(mesh);
    vacuum_interpolator.extract_solution(&laplace_solver);
    end_msg(t0);

    check_return(fields.check_limits(vacuum_interpolator.nodes.get_solutions()), "Field enhancement is out of limits!");

    vacuum_interpolator.nodes.write("out/result_E_phi.xyz");
    vacuum_interpolator.lintets.write("out/result_E_phi_linear.vtk");
    vacuum_interpolator.quadtets.write("out/result_E_phi_quad.vtk");

    laplace_solver.write("out/laplace.vtk");

    return fail;
}

} /* namespace femocs */
