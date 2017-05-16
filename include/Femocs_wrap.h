/* Are we compiling this with a C++ compiler? */
#ifdef __cplusplus 
extern "C" {
    class femocs::Femocs;
    typedef femocs::Femocs FEMOCS;

/* From the C side, we use an opaque pointer. */
#else
    typedef struct FEMOCS FEMOCS;
#endif

/* Constructor */
FEMOCS* create_femocs(const char* s);

/* Destructor */
void delete_femocs(FEMOCS* femocs);

void femocs_run(FEMOCS* femocs, int* retval, double E_field, const char* message);

void femocs_import_file(FEMOCS* femocs, int* retval, const char* s);

void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms, double* coordinates, double* box, int* nborlist);

void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms, double* x, double* y, double* z, int* types);

void femocs_export_atom_types(FEMOCS* femocs, int* retval, int n_atoms, int* types);

void femocs_export_elfield(FEMOCS* femocs, int* retval, int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

void femocs_export_temperature(FEMOCS* femocs, int* retval, int n_atoms, double* T);

void femocs_export_charge_and_force(FEMOCS* femocs, int* retval, int n_atoms, double* xq);

void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points, double* x, double* y,
        double* z, double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points, double* x, double* y, double* z, double* phi, int* flag);

void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg);

void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg);

#ifdef __cplusplus
}
#endif
