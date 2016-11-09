/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

#include <iostream>
#include <string>
#include <vector>

using namespace std;

/** If FILEWRITEMODE then file writes are operating normally, otherwise they return immediately. */
#define FILEWRITEMODE true

/** If ASSERTMODE then the asserts are operating. Disabling it makes the code to run faster. */
#define ASSERTMODE true

/** If VERBOSEMODE then the information about the code execution progress is printed to console. */
#define VERBOSEMODE true
              
/** Types of regions used in the simulation */
struct Types {
    const int NONE = 0;      //!< type of atom with unknown position
    const int BULK = 1;      //!< type of bulk material
    const int SURFACE = 2;   //!< type of open material surface
    const int VACUUM = 3;    //!< type of vacuum
    const int VACANCY = 3;   //!< type of vacancies
    const int PERIMETER = 4; //!< type of the rim/outer edge of surface
    const int FIXED = -1;    //!< type of fixed atoms
    const int XMIN = 5;      //!< type of atom on negative x-face of simulation cell
    const int YMIN = 6;      //!< type of atom on negative y-face of simulation cell
    const int ZMIN = 7;      //!< type of atom on negative z-face of simulation cell
    const int XMAX = 10;     //!< type of atom on positive x-face of simulation cell
    const int YMAX = 9;      //!< type of atom on positive y-face of simulation cell
    const int ZMAX = 8;      //!< type of atom on positive z-face of simulation cell
};

// Small hack to define Types only once
#ifdef MAINFILE
    Types TYPES;
#else
    extern Types TYPES;
#endif

/** Definition to handle cases where vital operation does not complete normally */
#define check_success(success, message) if (!(success)) { __success_fails(__FILE__, __LINE__, message); return; }

// Definitions for development in debug mode
#if ASSERTMODE
    /** Definition to give informative error if the requirement is not met */
    #define require(condition, message) if (!(condition)) __requirement_fails(__FILE__, __LINE__, message)

// In release(-like) versions nothing happens
#else
    #define require(condition, message) {}
#endif // ASSERTMODE

// Definitions for development in verbose mode
#if VERBOSEMODE
    /** Definition to give warning if the expectation is not met */
    #define expect(condition, message) if (!(condition)) __expectation_fails(__FILE__, __LINE__, message)
    
    /** Definition to print progress messages and to find the start time of code execution */
    #define start_msg(t0, message) t0 = __start_msg(message)

    /** Definition to print the execution time of code */
    #define end_msg(t0) __end_msg(t0);

// In release(-like) versions nothing happens
#else
    #define expect(condition, message) {}
    #define start_msg(t0, message) {}
    #define end_msg(t0) {}
#endif // VERBOSEMODE

/** Return mask of indices that are equal to the scalar */
const vector<bool> vector_equal(const vector<int> *v, const int s);

/** Return mask of indices that are not equal to the scalar */
const vector<bool> vector_not(const vector<int> *v, const int s);

/** Return mask of indices that are greater than the scalar */
const vector<bool> vector_greater(const vector<double> *v, const double s);

/** Return mask of indices that are greater or equal than the scalar */
const vector<bool> vector_greater_equal(const vector<double> *v, const double s);

/** Return mask of indices that are less than the scalar */
const vector<bool> vector_less(const vector<double> *v, const double s);

/** Return mask of indices that are less or equal than the scalar */
const vector<bool> vector_less_equal(const vector<double> *v, const double s);

/** Return sorting indexes for one vector */
vector<size_t> get_sort_indices(const vector<int> &v, const string& direction = "up");

/** Sum of the elements in vector */
const int vector_sum(const vector<bool> &v);
const int vector_sum(const vector<int> &v);
const double vector_sum(const vector<double> &v);

/** Determine whether the value is close to one of the boundary values or not */
const bool on_boundary(const double val, const double boundary1, const double boundary2, const double eps);

/** Determine whether the value is close to the boundary value or not */
const bool on_boundary(const double val, const double boundary, const double eps);

/** Extract file type from file name */
const string get_file_type(const string& file_name);

/** Throw an informative error if requirement fails */
void __requirement_fails(const char *file, int line, string message);

/** Throw an informative warning if expectation fails */
void __expectation_fails(const char *file, int line, string message);

/** Throw a informative message if some process fails */
void __success_fails(const char *file, int line, string message);

const double __start_msg(const char* message);

const void __end_msg(const double t0);

#endif /* MACROS_H_ */
