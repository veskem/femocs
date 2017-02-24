/*
 * Macros.h
 *
 *  Created on: 4.4.2016
 *      Author: veske
 */

#ifndef MACROS_H_
#define MACROS_H_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

/** If ASSERTMODE then the asserts are operating.
 * It must be controlled on compile time to enable deeper code optimisation. */
#define ASSERTMODE true

/** Types of regions used in the simulation */
struct Types {
    const int NONE = 0;      ///< type of atom with unknown position
    const int BULK = 1;      ///< type of bulk material
    const int SURFACE = 2;   ///< type of open material surface
    const int VACUUM = 3;    ///< type of vacuum
    const int VACANCY = 3;   ///< type of vacancies
    const int PERIMETER = 4; ///< type of the rim/outer edge of surface
    const int FIXED = -1;    ///< type of fixed atoms
    const int CLUSTER = -2;  ///< type of a cluster
    const int XMIN = 5;      ///< type of atom on negative x-face of simulation cell
    const int YMIN = 6;      ///< type of atom on negative y-face of simulation cell
    const int ZMIN = 7;      ///< type of atom on negative z-face of simulation cell
    const int XMAX = 10;     ///< type of atom on positive x-face of simulation cell
    const int YMAX = 9;      ///< type of atom on positive y-face of simulation cell
    const int ZMAX = 8;      ///< type of atom on positive z-face of simulation cell

    const int TETNODE = 1;      ///< node on the vertex of tetrahedron
    const int EDGECENTROID = 2; ///< node on the centroid of line
    const int FACECENTROID = 3; ///< node on the centroid of triangular face
    const int TETCENTROID = 4;  ///< node on the centroid of tetrahedron
};

/** Flags to control the output behaviour of the code */
struct Modes {
    bool QUIET = false;      ///< If QUIET no information about the code execution progress is printed to console.
    bool VERBOSE = true;       ///< If DEBUG all the information about the code execution progress is printed to console.
    bool WRITEFILE = true;   ///< If WRITEFILE then file writers operate normally, otherwise they return immediately.
    bool PERIODIC = true;    ///< Imported atoms have periodic boundaries in x- & y-direction
};

// Small hack to define Types and Modes only once
#ifdef MAINFILE
    Types TYPES;
    Modes MODES;
#else
    extern Types TYPES;
    extern Modes MODES;
#endif

// Asserts for catching errors in development mode
#if ASSERTMODE
    /** Definition to give informative error if the requirement is not met */
    #define require(condition, message) if (!(condition)) requirement_fails(__FILE__, __LINE__, message)

    /** Definition to give informative warning if the expectation is not met */
    #define expect(condition, message)  if (!(condition)) expectation_fails(__FILE__, __LINE__, message)

// In release(-like) versions nothing happens
#else
    #define require(condition, message) {}
    #define expect(condition, message) {}
#endif // ASSERTMODE

/** Definition to handle cases where operation does not complete normally */
#define check_message(condition, message) if (condition) { write_log(message); if (!MODES.QUIET) cout << "\nFEMOCS: " << message << endl; return 1; }

/** Definition to print progress messages and to find the start time of code execution */
#define start_msg(t0, message) { write_log(message); if (MODES.VERBOSE) t0 = __start_msg(message); }

/** Definition to print the execution time of code */
#define end_msg(t0) if (MODES.VERBOSE) __end_msg(t0)

/** Return mask of indices that are equal to the scalar */
vector<bool> vector_equal(const vector<int> *v, const int s);

/** Return mask of indices that are not equal to the scalar */
vector<bool> vector_not(const vector<int> *v, const int s);

/** Return mask of indices that are greater than the scalar */
vector<bool> vector_greater(const vector<double> *v, const double s);

/** Return mask of indices that are greater or equal than the scalar */
vector<bool> vector_greater_equal(const vector<double> *v, const double s);
vector<bool> vector_greater_equal(const vector<int> *v, const int s);

/** Return mask of indices that are less than the scalar */
vector<bool> vector_less(const vector<double> *v, const double s);

/** Return mask of indices that are less or equal than the scalar */
vector<bool> vector_less_equal(const vector<double> *v, const double s);

/** Return sorting indexes for one vector */
vector<int> get_sort_indices(const vector<int> &v, const string& direction="up");

/** Sum of the elements in vector */
int vector_sum(const vector<bool> &v);
int vector_sum(const vector<int> &v);
double vector_sum(const vector<double> &v);

/** Determine whether the value is close to one of the boundary values or not */
bool on_boundary(const double val, const double boundary1, const double boundary2, const double eps);

/** Determine whether the value is close to the boundary value or not */
bool on_boundary(const double val, const double boundary, const double eps);

/** Extract file type from file name */
string get_file_type(const string& file_name);

/** Throw an informative error if requirement fails */
void requirement_fails(const char *file, int line, string message);

/** Throw an informative warning if expectation fails */
void expectation_fails(const char *file, int line, string message);

void write_message(const string& message);

/** Append line to log-file */
void write_log(const string& message);

double __start_msg(const char* message);

void __end_msg(const double t0);

#endif /* MACROS_H_ */
