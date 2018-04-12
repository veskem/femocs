/*
 * Constants.h
 *
 *  Created on: 9.2.2018
 *      Author: veske
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>

using namespace std;
namespace femocs {

enum BoundaryId {
    copper_surface = 1,
    vacuum_top = 2,
    copper_bottom = 3,
    vacuum_sides = 4,
    copper_sides = 5
};

/** Types of regions used in the simulation */
struct Types {
    const int NONE = 0;      ///< type of atom with unknown position
    const int BULK = 1;      ///< type of bulk material
    const int SURFACE = 2;   ///< type of open material surface
    const int VACUUM = 3;    ///< type of vacuum
    const int VACANCY = 3;   ///< type of vacancies
    const int PERIMETER = 4; ///< type of the rim/outer edge of surface
    const int TOP = 5;   ///< type of vacancies
    const int BOTTOM = 6; ///< type of the rim/outer edge of surface
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
    bool MUTE = false;       ///< If QUIET no information about the code execution progress is printed to console.
    bool VERBOSE = true;     ///< If VERBOSE all the information about the code execution progress is printed to console.
    bool WRITEFILE = true;   ///< If WRITEFILE then file writers operate normally, otherwise they return immediately.
    bool WRITELOG = true;    ///< If WRITELOG then writing log file is enabled
    bool PERIODIC = true;    ///< Imported atoms have periodic boundaries in x- & y-direction
};

struct Globals {
    double TIME = 0;         ///< Simulation time in fs
};

// Small hack to define structs only once
#ifdef MAINFILE
    Types TYPES;
    Modes MODES;
    Labels LABELS;
    Globals GLOBALS;
#else
    extern Types TYPES;
    extern Modes MODES;
    extern Labels LABELS;
    extern Globals GLOBALS;
#endif

} /* namespace femocs */

#endif /* GLOBALS_H_ */
