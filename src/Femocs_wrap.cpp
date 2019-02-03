
#define MAINFILE

#include "Femocs.h"
#include "Femocs_wrap.h"
#include <iostream>

FEMOCS* create_femocs(const char* s) {
    return new femocs::Femocs(string(s));
}

void delete_femocs(FEMOCS* femocs) {
    delete femocs;
}

void femocs_run(FEMOCS* femocs, int* retval, int timestep, double time) {
    retval[0] = femocs->run(timestep, time);
}

void femocs_import_file(FEMOCS* femocs, int* retval, const char* file_name) {
    retval[0] = femocs->import_atoms(string(file_name));
}

void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms,
        const double* data, const double* box, const int* nborlist)
{
    retval[0] = femocs->import_parcas(n_atoms, data, box, nborlist);
}

void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms,
        const double* x, const double* y, const double* z, const int* types)
{
    retval[0] = femocs->import_atoms(n_atoms, x, y, z, types);
}

void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag)
{
    retval[0] = femocs->interpolate_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);
}

void femocs_interpolate_surface_elfield(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag)
{
    retval[0] = femocs->interpolate_surface_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);
}

void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z, double* phi, int* flag)
{
    retval[0] = femocs->interpolate_phi(n_points, x, y, z, phi, flag);
}

void femocs_export_data(FEMOCS* femocs, int* retval, double* data,
        int n_points, const char* data_type)
{
    retval[0] = femocs->export_data(data, n_points, data_type);
}

void femocs_export_int(FEMOCS* femocs, int* retval, int* data,
        int n_points, const char* data_type)
{
    retval[0] = femocs->export_data(data, n_points, data_type);
}

void femocs_interpolate(FEMOCS* femocs, int* retval, double* data, int* flag,
        int n_points, const char* data_type, int near_surface,
        const double* x, const double* y, const double* z)
{
    retval[0] = femocs->interpolate(data, flag, n_points, data_type, near_surface, x, y, z);
}

void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

void femocs_parse_string(FEMOCS* femocs, int* retval, const char* command, char* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}
