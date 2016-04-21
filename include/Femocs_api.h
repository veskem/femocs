
#ifndef FEMOCS_API_H_
#define FEMOCS_API_H_

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

// The const qualificators maps from the member function to pointers to the class instances.
const void femocs_run(FEMOCS* femocs, double d);

// Standalone function to call Femocs
const void femocs_speaker(const char* s);

#ifdef __cplusplus
}
#endif

#endif /* FEMOCS_API_H_ */
