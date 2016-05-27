/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

/** In case of debugmode == true more information is printed out and all the asserts are operating.
 * Disabling the debug mode makes code to run faster and removes the debug messages */
#define DEBUGMODE true

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
//namespace femocs {

#define require(condition, message) \
        if (!(condition))              \
          __requirement_fails(__FILE__, __LINE__, message)

// Definitions for development in debugging mode
#if DEBUGMODE
#define expect(condition, message) \
        if (!(condition))              \
          __expectation_fails(__FILE__, __LINE__, message)

/** Definition to print progress messages and to find the start time of code execution */
#define start_msg(t0, message) t0 = __start_msg(message)

/** Definition to print the execution time of code */
#define end_msg(t0) __end_msg(t0);

// In release(-like) versions nothing happens
#else
#define expect(condition, message) {}
#define start_msg(t0, message) {}
#define end_msg(t0) {}
#endif // DEBUGMODE

/** Sum of the elements in vector */
const inline int vector_sum(vector<bool> v) { return accumulate(v.begin(), v.end(), 0); }
const inline int vector_sum(vector<int> v) { return accumulate(v.begin(), v.end(), 0); }
const inline double vector_sum(vector<double> v) { return accumulate(v.begin(), v.end(), 0); }

/** Return mask of indices that doesn't contain the entry */
const vector<bool> vector_not(const vector<int> *v, const int entry);

/** Return mask of indices that do contain the entry */
const vector<bool> vector_equal(vector<int> *v, const int entry);

/** Determine whether the value is close to one of the boundary values or not */
const inline bool on_boundary(const double val, const double boundary1, const double boundary2, const double eps) {
    return (fabs(val - boundary1) <= eps) || (fabs(val - boundary2) <= eps);
}

/** Determine whether the value is close to the boundary value or not */
const inline bool on_boundary(const double val, const double boundary, const double eps) {
    return fabs(val - boundary) <= eps;
}

/** Extract file type from file name */
const inline string get_file_type(const string file_name) {
    const int start = file_name.find_last_of('.') + 1;
    const int end = file_name.size();
    return file_name.substr(start, end);
}

/** Throw an informative exception if expectation fails */
void __requirement_fails(const char *file, int line, string message);

/** Throw an warning if expectation fails */
void __expectation_fails(const char *file, int line, string message);

const double __start_msg(const char* message);

const void __end_msg(const double t0);

//} /* namespace femocs */

#endif /* MACROS_H_ */
