/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

/** In case of debugmode = true more information is printed out and all the asserts are operating.
 * Disabling the debug mode makes code to run faster and removes the debug messages */
#define DEBUGMODE true

#define STANDALONEMODE false
#define HELMODMODE true
#define KIMOCSMODE false

#include <string>
#include <vector>

using namespace std;
namespace femocs {

// Definitions for development and debugging mode
#if DEBUGMODE
#define assert(condition, message) \
        if (!(condition))              \
          requirement_fails(__FILE__, __LINE__, message)

/** Definition to print progress messages and to find the start time of code execution */
#define start_msg(t0, message) t0 = __start_msg(message)

/** Definition to print the execution time of code */
#define end_msg(t0) __end_msg(t0);

// In release(-like) versions nothing happens
#else
#define assert(condition, message) {}
#define start_msg(t0, message) {}
#define end_msg(t0) {}
#endif // DEBUGMODE

/** Throw an informative exception if expectation fails */
void requirement_fails(const char *file, int line, char* message);

const double __start_msg(const char* message);

const void __end_msg(const double t0);

} /* namespace femocs */

#endif /* MACROS_H_ */
