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
//namespace femocs {

/** If ASSERTMODE then the asserts are operating.
 * It must be controlled on compile time to enable deeper code optimisation. */
#define ASSERTMODE true

/** Throw an informative error if requirement fails */
void requirement_fails(const char *file, int line, string message);

/** Throw an informative warning if expectation fails */
void expectation_fails(const char *file, int line, string message);

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
#define check_return(condition, message) if (condition) { write_silent_msg(message); return 1; }

/** Return mask of indices that are equal to the scalar */
vector<bool> vector_equal(const vector<int> *v, const int s);

/** Return mask of indices that are not equal to the scalar */
vector<bool> vector_not(const vector<int> *v, const int s);

/** Return mask of indices that are greater than the scalar */
vector<bool> vector_greater(const vector<double> *v, const double s);
vector<bool> vector_greater(const vector<int> *v, const int s);

/** Return mask of indices that are greater or equal than the scalar */
vector<bool> vector_greater_equal(const vector<double> *v, const double s);
vector<bool> vector_greater_equal(const vector<int> *v, const int s);

/** Return mask of indices that are less than the scalar */
vector<bool> vector_less(const vector<double> *v, const double s);
vector<bool> vector_less(const vector<int> *v, const int s);

/** Return mask of indices that are less or equal than the scalar */
vector<bool> vector_less_equal(const vector<double> *v, const double s);
vector<bool> vector_less_equal(const vector<int> *v, const int s);

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

/** Write message to log file and if in silent or verbose mode, also to console */
void write_silent_msg(const string& message);

/** Write message to log file and if in verbose mode, also to console */
void write_verbose_msg(const string& message);

/** Periodic image of a point with respect to box xmin xmax */
double periodic_image(double p, double max, double min);

/** Append line to log-file */
void write_log(const string& message);

/** Clear the contents of log-file */
void clear_log();

/** Print progress messages and find the start time of code execution */
void start_msg(double& t0, const string& message);

/** Print the execution time of the code */
void end_msg(const double t0);

//} /* namespace femocs */
#endif /* MACROS_H_ */
