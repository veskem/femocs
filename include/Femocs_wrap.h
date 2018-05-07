/* Are we compiling this with a C++ compiler? */
#ifdef __cplusplus 
extern "C" {
    typedef femocs::Femocs FEMOCS;

/* From the C side, we use an opaque pointer. */
#else
    typedef struct FEMOCS FEMOCS;
#endif

/* Constructor */
FEMOCS* create_femocs(const char* s);

/* Destructor */
void delete_femocs(FEMOCS* femocs);

void femocs_run(FEMOCS* femocs, int* retval, int timestep);

void femocs_import_file(FEMOCS* femocs, int* retval, const char* file_name);

void femocs_import_parcas(FEMOCS* femocs, int* retval, int n_atoms,
        const double* coordinates, const double* box, const int* nborlist);

void femocs_import_atoms(FEMOCS* femocs, int* retval, int n_atoms,
        const double* x, const double* y, const double* z, const int* types);

void femocs_interpolate_elfield(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

void femocs_interpolate_surface_elfield(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z,
        double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

void femocs_interpolate_phi(FEMOCS* femocs, int* retval, int n_points,
        const double* x, const double* y, const double* z,
        double* phi, int* flag);

void femocs_export_data(FEMOCS* femocs, int* retval, double* data, int n_points, const char* data_type);

void femocs_interpolate(FEMOCS* femocs, int* retval, double* data, int* flag, int n_points, const char* data_type,
        int near_surface, const double* x, const double* y, const double* z);

void femocs_parse_int(FEMOCS* femocs, int* retval, const char* command, int* arg);

void femocs_parse_double(FEMOCS* femocs, int* retval, const char* command, double* arg);

void femocs_parse_string(FEMOCS* femocs, int* retval, const char* command, char* arg);

#ifdef __cplusplus
}
#endif
