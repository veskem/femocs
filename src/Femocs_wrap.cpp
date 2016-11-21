
#define MAINFILE

#include "Femocs.h"
#include <iostream>

using namespace std;

// ====================== Header ========================

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
    class femocs::Femocs;
    typedef femocs::Femocs FEMOCS;
#else
    // From the C side, we use an opaque pointer.
    typedef struct FEMOCS FEMOCS;
#endif

// Constructor
FEMOCS* create_femocs(const char* s);

// Destructor
void delete_femocs(FEMOCS* femocs);

const void femocs_run(FEMOCS* femocs, int* retval, double E_field, const char* message);

const void femocs_import_file(FEMOCS* femocs, int* retval, const char* s);

const void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms, double* coordinates, double* box, int* nborlist);

const void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms, double* x, double* y, double* z, int* types);

const void femocs_export_elfield(FEMOCS* femocs, int* retval, int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

const void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

const void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z, double* phi, int* flag);

const void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg);

const void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg);

#ifdef __cplusplus
}
#endif

// =================== Implementation =======================

FEMOCS* create_femocs(const char* s){
    return new femocs::Femocs(string(s));
}


void delete_femocs(FEMOCS* femocs){
    delete femocs;
}


const void femocs_run(FEMOCS* femocs, int* retval, double E_field, const char* message){
    retval[0] = femocs->run(E_field, string(message));
}


const void femocs_import_file(FEMOCS* femocs, int* retval, const char* s) {
    retval[0] = femocs->import_atoms(string(s));
}


const void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms, double* coordinates, double* box, int* nborlist) {
    retval[0] = femocs->import_atoms(n_atoms, coordinates, box, nborlist);
}


const void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms, double* x, double* y, double* z, int* types){
    retval[0] = femocs->import_atoms(n_atoms, x, y, z, types);
}


const void femocs_export_elfield(FEMOCS* femocs, int* retval, int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm){
    retval[0] = femocs->export_elfield(n_atoms, Ex, Ey, Ez, Enorm);
}


const void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag) {
    retval[0] = femocs->interpolate_elfield(n_points, x, y, z, Ex, Ey, Ez, Enorm, flag);
}


const void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z, double* phi, int* flag) {
    retval[0] = femocs->interpolate_phi(n_points, x, y, z, phi, flag);
}


const void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}


const void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg) {
    retval[0] = femocs->parse_command(string(command), arg);
}

