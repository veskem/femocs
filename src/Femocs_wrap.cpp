
#include "Femocs.h"
#include <iostream>

using namespace std;

// ====================== Header ========================

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
    class Femocs;
    typedef Femocs FEMOCS;
#else
    // From the C side, we use an opaque pointer.
    typedef struct FEMOCS FEMOCS;
#endif

// Constructor
FEMOCS* create_femocs(const char* s);

// Destructor
void delete_femocs(FEMOCS* femocs);

const void femocs_run(FEMOCS* femocs, double E_field, double*** phi);

const void femocs_import_atoms(FEMOCS* femocs, int n_atoms, double* x, double* y, double* z, int* types);

// Standalone function to call Femocs
const void femocs_speaker(const char* s);

#ifdef __cplusplus
}
#endif

// =================== Implementation =======================

FEMOCS* create_femocs(const char* s){
    return new Femocs(string(s));
}

void delete_femocs(FEMOCS* femocs){
    delete femocs;
}

const void femocs_run(FEMOCS* femocs, double E_field, double*** phi){
    femocs->run(E_field, phi);
}

const void femocs_import_atoms(FEMOCS* femocs, int n_atoms, double* x, double* y, double* z, int* types){
    femocs->import_atoms(n_atoms, x, y, z, types);
}

const void femocs_speaker(const char* s) {
    femocs_speaker(string(s));
}
