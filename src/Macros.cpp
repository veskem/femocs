/*
 * Macros.cpp
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#include "Macros.h"

using namespace std;
namespace femocs {

void requirement_fails(const char *file, int line, string message) {
    string exc = "Exception:\nfile = " + string(file) + "\nline = " + d2s(line)
            + "\nmessage = " + message + "\n";
    throw runtime_error(exc);
}

} /* namespace femocs */
