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
    const int EVAPORATED= -3;///< type of a evaporated atom
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

    /// Vtk cell types
    static struct Vtk_Types {
        // Linear cells
        static constexpr int EMPTY_CELL       = 0;
        static constexpr int VERTEX           = 1;
        static constexpr int POLY_VERTEX      = 2;
        static constexpr int LINE             = 3;
        static constexpr int POLY_LINE        = 4;
        static constexpr int TRIANGLE         = 5;
        static constexpr int TRIANGLE_STRIP   = 6;
        static constexpr int POLYGON          = 7;
        static constexpr int PIXEL            = 8;
        static constexpr int QUADRANGLE       = 9;
        static constexpr int TETRAHEDRON      = 10;
        static constexpr int VOXEL            = 11;
        static constexpr int HEXAHEDRON       = 12;
        static constexpr int WEDGE            = 13;
        static constexpr int PYRAMID          = 14;
        static constexpr int PENTAGONAL_PRISM = 15;
        static constexpr int HEXAGONAL_PRISM  = 16;
        // Quadratic, isoparametric cells
        static constexpr int QUADRATIC_EDGE                   = 21;
        static constexpr int QUADRATIC_TRIANGLE               = 22;
        static constexpr int QUADRATIC_QUADRANGLE             = 23;
        static constexpr int QUADRATIC_POLYGON                = 36;
        static constexpr int QUADRATIC_TETRAHEDRON            = 24;
        static constexpr int QUADRATIC_HEXAHEDRON             = 25;
        static constexpr int QUADRATIC_WEDGE                  = 26;
        static constexpr int QUADRATIC_PYRAMID                = 27;
        static constexpr int BIQUADRATIC_QUADRANGLE           = 28;
        static constexpr int TRIQUADRATIC_HEXAHEDRON          = 29;
        static constexpr int QUADRATIC_LINEAR_QUADRANGLE      = 30;
        static constexpr int QUADRATIC_LINEAR_WEDGE           = 31;
        static constexpr int BIQUADRATIC_QUADRATIC_WEDGE      = 32;
        static constexpr int BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33;
        static constexpr int BIQUADRATIC_TRIANGLE             = 34;
        // Polyhedron cell (consisting of polygonal faces)
        static constexpr int POLYHEDRON = 42;
    } VTK;
};

/** Labels of the exported/interpolated data.
 * Char component is used while calling Femocs externally, while string component is used internally.
 * Note that labels for vectors are capitalised, while lower case designates scalars.
 * While adding new labels, make sure not to use already taken char. */
struct Labels {
    const pair<char, string> vec = {'X', "vec"};
    const pair<char, string> vec_norm = {'y', "vec_norm"};
    const pair<char, string> scalar = {'z', "scalar"};
    const pair<char, string> elfield = {'E', "elfield"};
    const pair<char, string> elfield_norm = {'e', "elfield_norm"};
    const pair<char, string> potential = {'u', "potential"};
    const pair<char, string> temperature = {'t', "temperature"};
    const pair<char, string> rho = {'C', "rho"};
    const pair<char, string> rho_norm = {'c', "rho_norm"};
    const pair<char, string> pair_potential = {'p', "pair_potential"};
    const pair<char, string> pair_potential_sum = {'s', "pair_potential_sum"};
    const pair<char, string> force = {'F', "force"};
    const pair<char, string> force_norm = {'f', "force_norm"};
    const pair<char, string> charge = {'q', "charge"};
    const pair<char, string> velocity = {'V', "velocity"};
    const pair<char, string> velocity_norm = {'v', "velocity_norm"};

    int size() const { return (&velocity_norm)-(&vec) + 1; }  ///< number of values in Commands

    /** Convert given short char command into string label */
    string decode(const char c) {
        pair<char,string> const *ptr = &vec - 1;
        for (int i = 0; i < size(); ++i)
            if ((++ptr)->first == c)
                return ptr->second;
        return "";
    }

    /** Convert given string label into short char command */
    char encode(const string &s) {
        pair<char,string> const *ptr = &vec - 1;
        for (int i = 0; i < size(); ++i)
            if ((++ptr)->second == s)
                return ptr->first;
        return '\0';
    }
};

/** Flags to control the output behaviour of the code */
struct Modes {
    bool MUTE = false;      ///< If QUIET no information about the code execution progress is printed to console.
    bool VERBOSE = true;     ///< If VERBOSE all the information about the code execution progress is printed to console.
    bool WRITEFILE = true;   ///< If WRITEFILE then file writers operate normally, otherwise they return immediately.
    bool WRITELOG = true;    ///< If WRITELOG then writing log file is enabled
    bool PERIODIC = true;    ///< Imported atoms have periodic boundaries in x- & y-direction
};

// Small hack to define Types and Modes only once
#ifdef MAINFILE
    Types TYPES;
    Modes MODES;
    Labels LABELS;
#else
    extern Types TYPES;
    extern Modes MODES;
    extern Labels LABELS;
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

/** Throw an informative error if requirement fails */
void requirement_fails(const char *file, int line, string message);

/** Throw an informative warning if expectation fails */
void expectation_fails(const char *file, int line, string message);

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

#endif /* MACROS_H_ */
