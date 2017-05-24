
#define MAINFILE

#include "Femocs.h"
#include "Femocs_wrap.h"
#include <iostream>

FEMOCS* create_femocs(const char* s){
    return new femocs::Femocs(string(s));
}

void delete_femocs(FEMOCS* femocs){
    delete femocs;
}

void femocs_run(FEMOCS* femocs, int* retval, double E_field, const char* message){
    retval[0] = femocs->run(E_field, string(message));
}

void femocs_import_file(FEMOCS* femocs, int* retval, const char* s) {
    retval[0] = femocs->import_atoms(string(s));
}

void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms, double* coordinates, double* box, int* nborlist) {
    retval[0] = femocs->import_atoms(n_atoms, coordinates, box, nborlist);
}

void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms, double* x, double* y, double* z, int* types){
    retval[0] = femocs->import_atoms(n_atoms, x, y, z, types);
}

void femocs_export_atom_types(FEMOCS* femocs, int* retval, int n_atoms, int* types) {
    retval[0] = femocs->export_atom_types(n_atoms, types);
}

void femocs_export_elfield(FEMOCS* femocs, int* retval, int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm){
    retval[0] = femocs->export_elfield(n_atoms, Ex, Ey, Ez, Enorm);
}

void femocs_export_temperature(FEMOCS* femocs, int* retval, int n_atoms, double* T) {
    retval[0] = femocs->export_temperature(n_atoms, T);
}

void femocs_export_charge_and_force(FEMOCS* femocs, int* retval, int n_atoms, double* xq) {
    retval[0] = femocs->export_charge_and_force(n_atoms, xq);
}

void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    retval[0] = femocs->interpolate_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);
}

void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z, double* phi, int* flag) {
    retval[0] = femocs->interpolate_phi(n_points, x, y, z, phi, flag);
}

void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

void femocs_parse_boolean(FEMOCS* femocs, int* retval, const char* command, bool* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

void femocs_parse_string(FEMOCS* femocs, int* retval, const char* command, char* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}
