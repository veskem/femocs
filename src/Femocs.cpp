/*
 * Femocs.cpp
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#include "Femocs.h"
#include "Globals.h"
#include "ProjectRunaway.h"
#include "ProjectHeat.h"
#include "ProjectSpaceCharge.h"

#include <omp.h>

using namespace std;
namespace femocs {

Femocs::Femocs(const string &conf_file) : t0(0), reader(&conf.geometry) {
    static bool first_call = true;
    bool fail;

    // Read configuration parameters from configuration file
    conf.read_all(conf_file);

    // Pick the correct verbosity mode flags
    if      (conf.behaviour.verbosity == "mute")    { MODES.MUTE = true;  MODES.VERBOSE = false; }
    else if (conf.behaviour.verbosity == "silent")  { MODES.MUTE = false; MODES.VERBOSE = false; }
    else if (conf.behaviour.verbosity == "verbose") { MODES.MUTE = false; MODES.VERBOSE = true;  }

    // Pick correct flags for writing log file
    MODES.WRITELOG = conf.behaviour.n_write_log != 0;
    MODES.SHORTLOG = conf.behaviour.n_write_log < 0;

    // Clear the results from previous run
    if (first_call && conf.run.output_cleaner) fail = execute("rm -rf out");
    fail = execute("mkdir -p out");
    first_call = false;

    write_verbose_msg("======= Femocs started! =======");

    omp_set_num_threads(conf.behaviour.n_omp_threads);

    // pick the project to run
    if (conf.behaviour.project == "runaway")
        project = new ProjectRunaway(reader, conf);
    else if(conf.behaviour.project == "heat")
        project = new ProjectHeat(reader, conf);
    else if(conf.behaviour.project == "space_charge")
        project = new ProjectSpaceCharge(reader, conf);
    else {
        require(false, "Unimplemented project: " + conf.behaviour.project);
    }

    if (conf.path.restart_file != "")
        project->restart(conf.path.restart_file);
}

Femocs::~Femocs() {
    delete project;
    write_verbose_msg("======= Femocs finished! =======");
}

int Femocs::run(const int timestep, const double time) {
    return project->run(timestep, time);
}

void Femocs::perform_full_analysis(const int* nborlist) {
    string debug_msg = "Performing coordination";
    if (conf.run.rdf) debug_msg += ", rdf";
    if (conf.run.cluster_anal) debug_msg += ", cluster";
    start_msg(t0, debug_msg + " analysis");

    if (conf.run.rdf)
        reader.calc_rdf_coordinations(nborlist);
    else
        reader.calc_coordinations(nborlist);

    if (conf.run.cluster_anal)
        reader.calc_clusters(nborlist);

    end_msg(t0);
    write_verbose_msg(d2s(reader));

    start_msg(t0, "Extracting atom types");
    reader.extract_types();
    end_msg(t0);
}

void Femocs::perform_pseudo_analysis() {
    start_msg(t0, "Calculating coords from atom types");
    reader.calc_pseudo_coordinations();
    end_msg(t0);
}

int Femocs::import_atoms(const string& file_name, const int add_noise) {
    clear_log();
    string file_type, fname;

    if (file_name == "") fname = conf.path.infile;
    else fname = file_name;

    bool system_changed = true;
    if (fname == "generate") {
        start_msg(t0, "Generating nanotip");
        reader.generate_nanotip(conf.geometry.height, conf.geometry.radius, conf.geometry.latconst);
    } else {
        start_msg(t0, "Importing atoms");
        system_changed = reader.import_file(fname, add_noise);
    }
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + d2s(reader.size()) );

    file_type = get_file_type(fname);
    if (system_changed && file_type == "xyz")
        perform_full_analysis(NULL);
    else if (system_changed)
        perform_pseudo_analysis();
    else
        reader.extract_types();

    reader.write("atomreader.ckx");
    return 0;
}

int Femocs::import_parcas(const int n_atoms, const double* x0, const double* x1,
        const double* box, const int* nborlist)
{
    clear_log();

    start_msg(t0, "Importing coordinates & velocities");
    bool system_changed = reader.import_parcas(n_atoms, x0, x1, box, conf);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + d2s(reader.size()) );

    if (system_changed)
        perform_full_analysis(nborlist);
    else
        reader.extract_types();

    reader.write("atomreader.ckx");
    return 0;
}

int Femocs::import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types) {
    clear_log();
    conf.run.surface_cleaner = false; // disable the surface cleaner for atoms with known types

    start_msg(t0, "Importing atoms");
    bool system_changed = reader.import_atoms(n_atoms, x, y, z, types);
    end_msg(t0);
    write_verbose_msg( "#input atoms: " + d2s(reader.size()) );

    if (system_changed)
        perform_pseudo_analysis();

    reader.write("atomreader.ckx");
    return 0;
}

int Femocs::interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag)
{
    double *fields = new double[3*n_points];
    int retval = project->interpolate(fields, flag, n_points, LABELS.elfield, true, x, y, z);

    for (int i = 0; i < n_points; ++i) {
        int I = 3*i;
        Ex[i] = fields[I];
        Ey[i] = fields[I+1];
        Ez[i] = fields[I+2];
        Enorm[i] = sqrt(fields[I]*fields[I] + fields[I+1]*fields[I+1] + fields[I+2]*fields[I+2]);
    }

    delete[] fields;
    return retval;
}

int Femocs::interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag)
{
    double *fields = new double[3*n_points];
    int retval = project->interpolate(fields, flag, n_points, LABELS.elfield, false, x, y, z);

    for (int i = 0; i < n_points; ++i) {
        int I = 3*i;
        Ex[i] = fields[I];
        Ey[i] = fields[I+1];
        Ez[i] = fields[I+2];
        Enorm[i] = sqrt(fields[I]*fields[I] + fields[I+1]*fields[I+1] + fields[I+2]*fields[I+2]);
    }

    delete[] fields;
    return retval;
}

int Femocs::interpolate_phi(const int n_points, const double* x, const double* y, const double* z,
        double* phi, int* flag)
{
    return project->interpolate(phi, flag, n_points, LABELS.potential, false, x, y, z);
}

int Femocs::export_data(double* data, const int n_points, const string& data_type) {
    return project->export_data(data, n_points, data_type);
}

int Femocs::export_data(int* data, const int n_points, const string& data_type) {
    return project->export_data(data, n_points, data_type);
}

int Femocs::interpolate(double* data, int* flag,
        const int n_points, const string &data_type, const bool near_surface,
        const double* x, const double* y, const double* z)
{
    return project->interpolate(data, flag, n_points, data_type, near_surface, x, y, z);
}

int Femocs::parse_command(const string& command, int* arg) {
    return conf.read_command(command, arg[0]);
}

int Femocs::parse_command(const string& command, double* arg) {
    return conf.read_command(command, arg[0]);
}

int Femocs::parse_command(const string& command, string& arg) {
    return conf.read_command(command, arg);
}

int Femocs::parse_command(const string& command, char* arg) {
    string string_arg;
    bool fail = conf.read_command(command, string_arg);
    if (!fail) string_arg.copy(arg, string_arg.length());
    return fail;
}

} // namespace femocs
