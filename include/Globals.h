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

enum {
    n_coordinates    = 3,  ///< # coordinates
    n_nodes_per_edge = 2,  ///< # nodes on an edge
    n_nodes_per_tri  = 3,  ///< # nodes on a triangle
    n_nodes_per_quad = 4,  ///< # nodes on a quadrangle
    n_nodes_per_tet  = 4,  ///< # nodes on a tetrahedron
    n_nodes_per_hex  = 8,  ///< # nodes on a hexahedron
    n_edges_per_tri  = 3,  ///< # edges on a triangle
    n_edges_per_quad = 4,  ///< # edges on a quadrangle
    n_edges_per_tet  = 6,  ///< # edges on a tetrahedron
    n_edges_per_hex  = 12, ///< # edges on a hexahedron
    n_tris_per_tet   = 4,  ///< # triangles on a tetrahedron
    n_tets_per_tri   = 2,  ///< # tetrahedra connected to a triangle
    n_hexs_per_quad  = 2,  ///< # hexahedra connected to a quadrangle
    n_hexs_per_tet   = 4,  ///< # hexahedra connected to a tetrahedron
    n_quads_per_tri  = 3,  ///< # quadrangles connected to a triangle
    n_quads_per_hex  = 6   ///< # quadrangles connected to a hexahedron
};

/** Types of regions used in the simulation */
struct Types {
    const int NONE      = 0;  ///< type of atom with unknown position
    const int BULK      = 1;  ///< type of bulk material
    const int SURFACE   = 2;  ///< type of open material surface
    const int VACUUM    = 3;  ///< type of vacuum
    const int VACANCY   = 3;  ///< type of vacancies
    const int PERIMETER = 4;  ///< type of the rim/outer edge of surface
    const int FIXED     =-1;  ///< type of fixed atoms
    const int CLUSTER   =-2;  ///< type of a cluster
    const int EVAPORATED=-3;  ///< type of a evaporated atom
    const int XMIN      = 5;  ///< type of atom on negative x-face of simulation cell
    const int YMIN      = 6;  ///< type of atom on negative y-face of simulation cell
    const int ZMIN      = 7;  ///< type of atom on negative z-face of simulation cell
    const int XMAX      = 10; ///< type of atom on positive x-face of simulation cell
    const int YMAX      = 9;  ///< type of atom on positive y-face of simulation cell
    const int ZMAX      = 8;  ///< type of atom on positive z-face of simulation cell

    const int TETNODE      = 1; ///< node on the vertex of tetrahedron
    const int EDGECENTROID = 2; ///< node on the centroid of line
    const int FACECENTROID = 3; ///< node on the centroid of triangular face
    const int TETCENTROID  = 4; ///< node on the centroid of tetrahedron
};

/** IDs of mesh boundaries */
namespace BoundaryID {
enum {
    copper_surface = 2,
    vacuum_top     = 8,
    copper_bottom  = 7,
    vacuum_sides   = 4,
    copper_sides   = 4
};
}

/** Vtk cell types */
namespace VtkType {
enum {
    // Linear cells
    empty_cell       = 0,
    vertex           = 1,
    poly_vertex      = 2,
    line             = 3,
    poly_line        = 4,
    triangle         = 5,
    triangle_strip   = 6,
    polygon          = 7,
    pixel            = 8,
    quadrangle       = 9,
    tetrahedron      = 10,
    voxel            = 11,
    hexahedron       = 12,
    wedge            = 13,
    pyramid          = 14,
    pentagonal_prism = 15,
    hexagonal_prism  = 16,

    // Quadratic, isoparametric cells
    quadratic_edge                   = 21,
    quadratic_triangle               = 22,
    quadratic_quadrangle             = 23,
    quadratic_polygon                = 36,
    quadratic_tetrahedron            = 24,
    quadratic_hexahedron             = 25,
    quadratic_wedge                  = 26,
    quadratic_pyramid                = 27,
    biquadratic_quadrangle           = 28,
    triquadratic_hexahedron          = 29,
    quadratic_linear_quadrangle      = 30,
    quadratic_linear_wedge           = 31,
    biquadratic_quadratic_wedge      = 32,
    biquadratic_quadratic_hexahedron = 33,
    biquadratic_triangle             = 34,

    // Polyhedron cell (consisting of polygonal faces)
    polyhedron = 42
};
}

/** Gmesh cell types */
namespace GmshType {
enum {
    vertex                  = 15, ///< 1-node point
    line                    = 1,  ///< 2-node line
    triangle                = 2,  ///< 3-node triangle
    quadrangle              = 3,  ///< 4-node quadrangle
    tetrahedron             = 4,  ///< 4-node tetrahedron
    hexahedron              = 5,  ///< 8-node hexahedron
    quadratic_edge          = 8,  ///< 3-node second order line
    quadratic_triangle      = 9,  ///< 6-node second order triangle
    quadratic_quadrangle    = 16, ///< 8-node second order quadrangle
    quadratic_tetrahedron   = 11, ///< 10-node second order tetrahedron
    quadratic_hexahedron    = 17, ///< 20-node second order hexahedron
    biquadratic_quadrangle  = 10, ///< 9-node second order quadrangle
    triquadratic_hexahedron = 12  ///< 27-node second order hexahedron
};
}

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
    const string kin_energy = "kin_energy";
    const string pot_energy = "pot_energy";
    const string pair_potential = "pair_potential";
    const string parcas_force = "parcas_force";
    const string charge_force = "charge_force";
    const string force = "force";
    const string force_norm = "force_norm";
    const string charge = "charge";
    const string charge_density = "charge_density";
    const string parcas_velocity = "parcas_velocity";
    const string velocity = "velocity";
    const string velocity_norm = "velocity_norm";
    const string heat = "heat";
    const string area = "area";
    const string atom_type = "atom_type";
    const string nodes = "nodes";
    const string edges = "edges";
    const string triangles = "triangles";
    const string quadrangles = "quadrangles";
    const string tetrahedra = "tetrahedra";
    const string hexahedra = "hexahedra";
};

/** Flags to control the output behaviour of the code */
struct Modes {
    bool MUTE = false;       ///< If QUIET no information about the code execution progress is printed to console.
    bool VERBOSE = true;     ///< If VERBOSE all the information about the code execution progress is printed to console.
    bool WRITELOG = true;    ///< If WRITELOG then writing log file is enabled
    bool SHORTLOG = true;    ///< If SHORTLOG then only the last timestep is stored in log file
    bool PERIODIC = true;    ///< If PERIODIC then imported atoms have periodic boundaries in x- & y-direction
    double WRITE_PERIOD = 0; ///< Min time in fs between two consequent writes to the same file; 0 enables writing at every call and <0 turns writing off
    string OUT_FOLDER = "out";  ///< Path to the folder where files will be written
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
