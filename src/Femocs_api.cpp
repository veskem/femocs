
#include "Femocs_api.h"
#include "Femocs.h"

#include <iostream>

using namespace std;

FEMOCS* create_femocs(const char* s){
    return new Femocs(string(s));
}

void delete_femocs(FEMOCS* femocs){
    delete femocs;
}

const void femocs_run(FEMOCS* femocs, double d){
    femocs->run(d);
}

const void femocs_speaker(const char* s) {
    femocs_speaker(string(s));
}
