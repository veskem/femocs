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

enum MeshData {
    n_coordinates = 3,     ///< # coordinates
    n_nodes_per_edge = 2, ///< # nodes on an edge
    n_nodes_per_tri = 3,  ///< # nodes on a triangle
    n_nodes_per_quad = 4, ///< # nodes on a quadrangle
    n_nodes_per_tet = 4,  ///< # nodes on a tetrahedron
    n_nodes_per_hex = 8,  ///< # nodes on a hexahedron
    n_edges_per_tri = 3,  ///< # edges on a triangle
    n_edges_per_quad = 4, ///< # edges on a quadrangle
    n_edges_per_tet = 6,  ///< # edges on a tetrahedron
    n_edges_per_hex = 12, ///< # edges on a hexahedron
    n_tris_per_tet = 4,   ///< # triangles on a tetrahedron
    n_tets_per_tri = 2,   ///< # tetrahedra connected to a triangle
    n_hexs_per_quad = 2,  ///< # hexahedra connected to a quadrangle
    n_hexs_per_tet = 4,   ///< # hexahedra connected to a tetrahedron
    n_quads_per_tri = 3,  ///< # quadrangles connected to a triangle
    n_quads_per_hex = 6  ///< # quadrangles connected to a hexahedron
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

/** Labels of the calculated/exported/interpolated data. */
struct Labels {
    const string vec = "vec";
    const string vec_norm = "vec_norm";
    const string scalar = "scalar";
    const string elfield = "elfield";
    const string elfield_norm = "elfield_norm";
    const string potential = "potential";
    const string temperature = "temperature";
    const string rho = "rho";
    const string rho_norm = "rho_norm";
    const string pair_potential = "pair_potential";
    const string pair_potential_sum = "pair_potential_sum";
    const string parcas_force = "parcas_force";
    const string charge_force = "charge_force";
    const string force = "force";
    const string force_norm = "force_norm";
    const string charge = "charge";
    const string parcas_velocity = "parcas_velocity";
    const string velocity = "velocity";
    const string velocity_norm = "velocity_norm";
    const string heat = "heat";
    const string area = "area";
    const string atom_type = "atom_type";
};

/** Flags to control the output behaviour of the code */
struct Modes {
    bool MUTE = false;       ///< If QUIET no information about the code execution progress is printed to console.
    bool VERBOSE = true;     ///< If VERBOSE all the information about the code execution progress is printed to console.
    bool WRITEFILE = true;   ///< If WRITEFILE then file writers operate normally, otherwise they return immediately.
    bool WRITELOG = true;    ///< If WRITELOG then writing log file is enabled
    bool SHORTLOG = true;    ///< If SHORTLOG then only the last timestep is stored in log file
    bool PERIODIC = true;    ///< If PERIODIC then imported atoms have periodic boundaries in x- & y-direction
};

struct Globals {
    double TIME = 0;         ///< Simulation time in fs
    int TIMESTEP = 0;        ///< Simulation time step
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
