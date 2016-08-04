/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

/** If DEBUGMODE then the asserts are operating. Disabling it makes the code to run faster. */
#define DEBUGMODE true

/** If VERBOSE then the debug information about the code execution is printed out. */
#define VERBOSE true

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
//namespace femocs {

/** Types of regions used in the simulation */
struct Types {
    const int BULK = 1;    //!< type of bulk material
    const int SURFACE = 2; //!< type of open material surface
    const int VACANCY = 3; //!< type of vacancies
    const int VACUUM = 3;  //!< type of vacuum
    const int EDGE = 0;    //!< type of the rim/outer edge of surface
    const int FIXED = -1;  //!< type of fixed atoms
    const int XMIN = 4;    //!< type of atom on negative x-face of simulation cell
    const int YMIN = 5;    //!< type of atom on negative y-face of simulation cell
    const int ZMIN = 6;    //!< type of atom on negative z-face of simulation cell
    const int XMAX = 10;   //!< type of atom on positive x-face of simulation cell
    const int YMAX = 9;    //!< type of atom on positive y-face of simulation cell
    const int ZMAX = 8;    //!< type of atom on positive z-face of simulation cell
    const int NONE = 7;    //!< type of atom with unknown position
};

#ifdef MAINFILE
    Types TYPES;
#else
    extern Types TYPES;
#endif

#if DEBUGMODE
    #define require(condition, message) \
            if (!(condition)) __requirement_fails(__FILE__, __LINE__, message)

// In release(-like) versions nothing happens
#else
    #define require(condition, message) {}
#endif // DEBUGMODE

// Definitions for development in VERBOSE mode
#if VERBOSE
    /** Definition to print progress messages and to find the start time of code execution */
    #define start_msg(t0, message) t0 = __start_msg(message)

    /** Definition to print the execution time of code */
    #define end_msg(t0) __end_msg(t0);
    #define expect(condition, message) \
            if (!(condition))              \
              __expectation_fails(__FILE__, __LINE__, message)

// In release(-like) versions nothing happens
#else
    #define expect(condition, message) {}
    #define start_msg(t0, message) {}
    #define end_msg(t0) {}
#endif // VERBOSE

/** Sum of the elements in vector */
const inline int vector_sum(vector<bool> v) { return accumulate(v.begin(), v.end(), 0); }
const inline int vector_sum(vector<int> v) { return accumulate(v.begin(), v.end(), 0); }
const inline double vector_sum(vector<double> v) { return accumulate(v.begin(), v.end(), 0); }


/** Template to return mask of indices that satisfy the comparison condition with an entry */
template<typename T, typename Op>
const vector<bool> __vector_compare(const vector<T> *v, const T entry) {
    vector<bool> mask(v->size());
    Op op;
    for (int i = 0; i < v->size(); ++i)
        mask[i] = op((*v)[i], entry);

    return mask;
}

/** Return mask of indices that are equal to the scalar */
const vector<bool> vector_equal(const vector<int> *v, const int s);

/** Return mask of indices that are not equal to the scalar */
const vector<bool> vector_not(const vector<int> *v, const int s);

/** Return mask of indices that are greater than the scalar */
const vector<bool> vector_greater(const vector<double> *v, const double s);

/** Return mask of indices that are greater or equal than the scalar */
const vector<bool> vector_ge(const vector<double> *v, const double s);

/** Return mask of indices that are less than the scalar */
const vector<bool> vector_less(const vector<double> *v, const double s);

/** Return mask of indices that are less or equal than the scalar */
const vector<bool> vector_less_equal(const vector<double> *v, const double s);


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
