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
Femocs::Femocs(const string &conf_file) : skip_calculations(false), fail(false) {
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
    if (first_call && conf.clear_output) fail = system("rm -rf output");
    fail = system("mkdir -p output");

    start_msg(t0, "======= Femocs started! =======\n");

    // Initialise heating module
    start_msg(t0, "=== Reading physical quantities...");
    ch_solver1.set_physical_quantities(&phys_quantities);
    ch_solver2.set_physical_quantities(&phys_quantities);
    ch_solver  = &ch_solver1;
    prev_ch_solver = NULL;
    end_msg(t0);

    first_call = false;
}

// delete data and print bye-bye-message
Femocs::~Femocs() {
    start_msg(t0, "======= Femocs finished! =======\n");
}

// Workhorse function to generate FEM mesh and to solve differential equation(s)
int Femocs::run(const double elfield, const string &message) {
    static unsigned int timestep = 0;  // Counter to measure how many times Femocs has been called
    static bool prev_skip_calculations = true;  // Value of skip_calculations in last call
    double tstart;                     // Variable used to measure the code execution time
    stringstream stream; stream << fixed << setprecision(3);

    if (!prev_skip_calculations && MODES.WRITEFILE)
        MODES.WRITEFILE = false;

    if ((conf.n_writefile > 0) && (timestep % conf.n_writefile == 0))
        MODES.WRITEFILE = true;

    conf.message = to_string(timestep++);
    write_silent_msg("Running at timestep " + conf.message);
    conf.message = "_" + string( max(0.0, 5.0 - conf.message.length()), '0' ) + conf.message;

    prev_skip_calculations = skip_calculations;

    stream.str(""); stream << "Atoms haven't moved significantly, " << reader.rms_distance
        << " < " << conf.distance_tol << "! Field calculation will be skipped!";
    check_return(skip_calculations, stream.str());

    skip_calculations = true;
    conf.E0 = elfield; // long-range electric field
    conf.neumann = -10.0 * elfield;  // set minus gradient of solution to equal to E0; also convert V/Angstrom  to  V/nm
    tstart = omp_get_wtime();

    TetgenMesh bulk_mesh;   // FEM mesh in bulk material
    TetgenMesh vacuum_mesh; // FEM mesh in vacuum

    // Generate FEM mesh
    check_return( generate_meshes(bulk_mesh, vacuum_mesh), "Mesh generation failed!" );

    // Store parameters for comparing the results with analytical hemi-ellipsoid results
    vacuum_interpolator.set_analyt(conf.E0, conf.radius, dense_surf.sizes.zbox);
    fields.set_analyt(conf.E0, conf.radius, dense_surf.sizes.zbox);

    // Solve Laplace equation on vacuum mesh
    fch::Laplace<3> laplace_solver;
    check_return( solve_laplace(vacuum_mesh, laplace_solver), "Solving Laplace equation failed!");

    // Solve heat & continuity equation on bulk mesh
    check_return( solve_heat(bulk_mesh, laplace_solver), "Solving heat & continuity equation failed!");

    // Extract face charges
    if( extract_charge(vacuum_mesh) ) {
        MODES.WRITEFILE = true;
        bulk_mesh.hexahedra.write("output/hexmesh_bulk" + conf.message + ".vtk");
        vacuum_mesh.faces.write("output/surface_faces_clean.vtk");
        vacuum_interpolator.write("output/result_E_phi.vtk");
        vacuum_interpolator.write("output/result_E_phi.xyz");
        face_charges.write("output/face_charges.vtk");
        check_return(true, "Error calculating face charges!");
    }

    start_msg(t0, "=== Saving atom positions...");
    reader.save_current_run_points(conf.distance_tol);
    end_msg(t0);

    stream.str(""); stream << "Total execution time " << omp_get_wtime() - tstart;
    write_silent_msg(stream.str());
    skip_calculations = false;

//    write_slice("output/slice.xyz");

    return 0;
}

// Generate boundary nodes for mesh
int Femocs::generate_boundary_nodes(Media& bulk, Media& coarse_surf, Media& vacuum) {
    start_msg(t0, "=== Extracting surface...");
    dense_surf.extract(reader, TYPES.SURFACE);
    dense_surf.sort_atoms(3, "down");
    end_msg(t0);

    dense_surf.write("output/surface_dense.xyz");

    if (conf.surface_cleaner == "voronois") {
        start_msg(t0, "=== Cleaning surface with Voronoi cells...");
//        fail = dense_surf.voronoi_clean(areas, conf.radius, conf.latconst, conf.mesh_quality + "a10");
        fail = dense_surf.voronoi_clean(areas, conf.radius, conf.latconst, conf.mesh_quality);
        check_return(fail, "Making voronoi cells failed!");
        end_msg(t0);

        dense_surf.write("output/surface_dense_clean.xyz");
    }

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
    coarse_surf = extended_surf;
    coarse_surf += dense_surf;
    coarse_surf = coarse_surf.clean(coarseners);
    coarse_surf.smoothen(conf.radius, conf.surface_smooth_factor, 3.0*conf.coord_cutoff);
    end_msg(t0);

    coarse_surf.write("output/surface_coarse.xyz");

    start_msg(t0, "=== Generating bulk & vacuum...");
    coarse_surf.calc_statistics();  // calculate zmin and zmax for surface
    vacuum.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin + conf.box_height * coarse_surf.sizes.zbox);
    bulk.generate_simple(coarse_surf.sizes, coarse_surf.sizes.zmin - conf.bulk_height * conf.latconst);
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
    // r - reconstruct, n - output neighbour list, Q - quiet, q - mesh quality, a - element volume,
    // F - suppress output of faces and edges, B - suppress output of boundary info
    string command = "rnQFBq" + conf.mesh_quality;
    if (conf.element_volume != "") command += "a" + conf.element_volume;
    fail = big_mesh.generate(bulk, coarse_surf, vacuum, command);

    check_return(fail, "Triangulation failed!");
    end_msg(t0);

    start_msg(t0, "=== Making surface faces...");
    big_mesh.generate_appendices();
    end_msg(t0);

    big_mesh.faces.write("output/surface_faces.vtk");

    start_msg(t0, "=== Marking tetrahedral mesh...");
    fail = big_mesh.mark_mesh();
    big_mesh.nodes.write("output/tetmesh_nodes.xyz");
    big_mesh.elems.write("output/tetmesh_elems.vtk");
    check_return(fail, "Mesh marking failed!");
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
    vacuum_mesh.faces.clean_sides(reader.sizes);
    end_msg(t0);

    vacuum_mesh.faces.write("output/surface_faces_clean.vtk");

    if (conf.surface_cleaner == "faces") {
        start_msg(t0, "=== Cleaning surface with triangles...");
        dense_surf.faces_clean(vacuum_mesh, conf.surface_thichness);
        end_msg(t0);

        dense_surf.write("output/surface_dense_clean.xyz");
    }

    expect(bulk_mesh.nodes.size() > 0, "Zero nodes in bulk mesh!");
    expect(vacuum_mesh.nodes.size() > 0, "Zero nodes in vacuum mesh!");
    expect(bulk_mesh.hexahedra.size() > 0, "Zero elements in bulk mesh!");
    expect(vacuum_mesh.hexahedra.size() > 0, "Zero elements in vacuum mesh!");

    bulk_mesh.hexahedra.write("output/hexmesh_bulk" + conf.message + ".vtk");
    vacuum_mesh.hexahedra.write("output/hexmesh_vacuum" + conf.message + ".vtk");

    stringstream ss; ss << "Bulk:   " << bulk_mesh << "\n  Vacuum: " << vacuum_mesh;
    write_verbose_msg(ss.str());

    return 0;
}

// Solve Laplace equation
int Femocs::solve_laplace(const TetgenMesh& mesh, fch::Laplace<3>& solver) {
    start_msg(t0, "=== Importing mesh to Laplace solver...");
    fail = !solver.import_mesh_directly(mesh.nodes.export_dealii(), mesh.hexahedra.export_dealii());
    check_return(fail, "Importing mesh to Deal.II failed!");
    end_msg(t0);

    start_msg(t0, "=== Initializing Laplace solver...");
    solver.set_applied_efield(conf.neumann);
    solver.setup_system();
    solver.assemble_system();
    end_msg(t0);

    stringstream ss; ss << solver;
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Running Laplace solver...");
    solver.solve(conf.n_phi, conf.phi_error, true, conf.ssor_param);
    end_msg(t0);

    start_msg(t0, "=== Extracting E and phi...");
    fail = vacuum_interpolator.extract_solution(&solver, mesh);
    end_msg(t0);

    vacuum_interpolator.write("output/result_E_phi.xyz");
    vacuum_interpolator.write("output/result_E_phi.vtk");
    vacuum_interpolator.print_statistics();
    vacuum_interpolator.print_enhancement();

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
    check_return(fail, "Importing mesh to Deal.II failed!");
    ch_solver->setup_system();
    end_msg(t0);

    stringstream ss; ss << *(ch_solver);
    write_verbose_msg(ss.str());

    start_msg(t0, "=== Transfering elfield to J & T solver...\n");
    FieldReader fr(&vacuum_interpolator);
    fr.interpolate(ch_solver, conf.use_histclean * conf.coord_cutoff, true);
    end_msg(t0);

    fr.write("output/surface_field.xyz");

    start_msg(t0, "=== Running J & T solver...\n");
    double t_error = ch_solver->run_specific(conf.t_error, conf.n_newton, false, "", MODES.VERBOSE, 2, 400, true);
    end_msg(t0);

    ch_solver->output_results("output/result_J_T.vtk");

    check_return(t_error > conf.t_error, "Temperature didn't converge, err=" + to_string(t_error));

    start_msg(t0, "=== Extracting T and rho...");
    bulk_interpolator.extract_solution(ch_solver, mesh);
    end_msg(t0);

    bulk_interpolator.write("output/result_J_T.xyz");

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

// Calculate the charges on surface faces
int Femocs::extract_charge(const TetgenMesh& mesh) {
    start_msg(t0, "=== Calculating face charges...");
    face_charges.calc_charges(mesh, conf.E0);             // electric field in the middle of face in directly from solution
    end_msg(t0);

    const double tot_charge = conf.E0 * reader.sizes.xbox * reader.sizes.ybox * face_charges.eps0;
    face_charges.print_statistics(tot_charge);

    check_return(!face_charges.charge_conserved(tot_charge, conf.charge_tolerance),
            "Face charges are not conserved!");

    face_charges.clean(dense_surf.sizes, conf.latconst);
    face_charges.write("output/face_charges.vtk");

    return 0;
}

// Interpolate the solution on the x-z plane in the middle of simulation box
void Femocs::write_slice(const string& file_name) {
	int writefile_save = MODES.WRITEFILE;
	MODES.WRITEFILE = true;

    const int nx = 300;  // number of points in x-direction
    const int nz = 300;  // number of points in z-direction
	const double eps = 1e-5;

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
    fr.interpolate(medium, conf.use_histclean * conf.coord_cutoff);
    fr.write(file_name);

   	MODES.WRITEFILE = writefile_save;
}

// import atoms from file
int Femocs::import_atoms(const string& file_name) {
    clear_log();
    string file_type, fname;

    if (file_name == "") fname = conf.atom_file;
    else fname = file_name;

    file_type = get_file_type(fname);
    require(file_type == "ckx" || file_type == "xyz", "Unknown file type: " + file_type);

    start_msg(t0, "=== Importing atoms...");
    reader.import_file(fname);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

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
            reader.check_clusters(1);

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
    write_verbose_msg( "#input atoms: " + to_string(reader.size()) );

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
        reader.check_clusters(1);

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

        reader.write("output/atomreader.xyz");
    }

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

    if (skip_calculations)
        write_silent_msg("Using previous solution!");
    else {
        start_msg(t0, "=== Interpolating E and phi...");
        fields.interpolate(dense_surf, conf.use_histclean * conf.coord_cutoff, 0, false);
        end_msg(t0);

        fields.write("output/fields.movie");
        fields.print_statistics();
        fields.print_enhancement();
    }

    start_msg(t0, "=== Exporting electric field...");
    fields.export_solution(n_atoms, Ex, Ey, Ez, Enorm);
    end_msg(t0);

    return 0;
}

// calculate and export temperatures on imported atom coordinates
int Femocs::export_temperature(const int n_atoms, double* T) {
    if (n_atoms < 0 || !conf.heating) return 0;
    check_return(bulk_interpolator.size() == 0, "No temperature to export!");

    if (skip_calculations)
        write_silent_msg("Using previous solution!");
    else {
        start_msg(t0, "=== Interpolating T and rho...");
        temperatures.interpolate(reader);
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
    if (n_atoms < 0) return 0;
    check_return(fields.size() == 0 || face_charges.size() == 0, "No force to export!");

    if (skip_calculations)
        write_silent_msg("Using previous solution!");
    else {
        start_msg(t0, "=== Calculating charges and forces...");
        forces.calc_forces(fields, face_charges, conf.use_histclean*conf.coord_cutoff,
                conf.charge_smooth_factor, conf.force_factor);

        if (conf.surface_cleaner == "voronois")
            forces.recalc_forces(fields, areas, conf.force_factor);
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
    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator.size() == 0, "No solution to export!");

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
int Femocs::interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
        double* phi, int* flag) {
    if (n_points <= 0) return 0;
    check_return(vacuum_interpolator.size() == 0, "No solution to export!");

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
