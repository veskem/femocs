/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

#include <sstream>
#include <stdexcept>
#include <string>

/** In case of debugmode = true more information is printed out and all the asserts are operating.
 * Disabling the debug mode makes code to run faster and removes the debug messages */
#define DEBUGMODE true

#define STANDALONEMODE true
#define HELMODMODE false
#define KIMOCSMODE false

using namespace std;
namespace femocs {

// Definitions for development and debugging mode
#if DEBUGMODE
#define assert(condition, message) \
        if (!(condition))              \
          requirement_fails(__FILE__, __LINE__, message)

/** Definition to print progress messages and to find the start time of code execution */
#define start_msg(t0, message)  \
        cout << endl << message;    \
        cout.flush();               \
        t0 = omp_get_wtime()

/** Definition to print the execution time of code */
#define end_msg(t0) cout << ", time: " << omp_get_wtime() - t0 << endl

// In release(-like) versions nothing happens
#else
#define assert(condition, message) {}
#define start_msg(t0, message) {}
#define end_msg(t0) {}
#endif // DEBUGMODE

/** Template to convert data to string */
template<typename T>
inline string d2s(T data) {
    ostringstream o;
    if (!(o << data)) throw runtime_error("Bad conversion of data to string!");
    return o.str();
}

/** Throw an informative exception if expectation fails */
void requirement_fails(const char *file, int line, string message);

} /* namespace femocs */
#endif /* MACROS_H_ */
