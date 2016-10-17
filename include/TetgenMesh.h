/*
 * TetgenMesh.h
 *
 *  Created on: 3.10.2016
 *      Author: veske
 */

#ifndef TETGENMESH_H_
#define TETGENMESH_H_

#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "TetgenCells.h"
#include "Media.h"

using namespace std;
namespace femocs {

/**
 * Class to create and handle finite element mesh in Tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/ */
class TetgenMesh {
public:
    TetgenMesh();
    ~TetgenMesh();

    /** Function to generate simple mesh that consists of one element */
    const void generate_simple();

    /** Function to generate mesh from surface, bulk and vacuum atoms */
    const void generate_mesh(Bulk &bulk, Surface &surf, Vacuum &vacuum, const string& cmd);

    const void generate_mesh_appendices();

    const bool mark_mesh(const bool postprocess);

    const void separate_meshes(TetgenMesh &vacuum, TetgenMesh &bulk, const string &cmd);

    const void write_tetgen(const string file_name);

    const void recalc();
    const void recalc(const string& cmd);
    const void recalc(const string& cmd1, const string& cmd2);

    // Tetgen data structure
    tetgenio tetIOin;
    tetgenio tetIOout;

    // Objects holding operations for accessing cell data
    TetgenNodes nodes = TetgenNodes(&tetIOout, &tetIOin);
    TetgenEdges edges = TetgenEdges(&tetIOout);
    TetgenFaces faces = TetgenFaces(&tetIOout);
    TetgenElements elems = TetgenElements(&tetIOout, &tetIOin);

    const int n_coordinates = 3;
    const int n_nodes_per_edge = 2;
    const int n_nodes_per_face = 3;
    const int n_nodes_per_elem = 4;
    const int n_edges_per_node = 3;
    const int n_edges_per_face = 3;
    const int n_edges_per_elem = 6;
    const int n_faces_per_elem = 4;

private:
    const void mark_nodes();
    const void mark_edges();
    const void mark_faces();
    const void mark_elems();

    const void remark_perimeter_nodes();
    const void remark_elems(const int skip_type);
    const bool post_process_marking();

    const void generate_edges();

    /** Function to generate surface faces from already existing elements and surface nodes */
    const void generate_surf_faces();
};

/** Class to mark Mesh nodes with ray-triangle intersection technique,
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle */
class RaySurfaceIntersect2 {
public:
    /** Constructor of RaySurfaceIntersect  */
    RaySurfaceIntersect2(TetgenMesh* mesh);

    /** Function to find with Moller-Trumbore algorithm whether the ray and the surface intersect or not */
    const bool ray_intersects_surface(const Vec3 &origin, const Vec3 &direction);

    /** Function to precompute the data needed to execute the Moller-Trumbore algorithm */
    const void precompute_triangles(const Vec3 &direction);

private:
    /** Constants to specify the tolerances */
    const double epsilon = 1e-3;
    const double zero = -1.0 * epsilon;
    const double one = 1.0 + epsilon;

    /** Pointer to Mesh with nodes and surface faces */
    TetgenMesh* mesh;

    /** Data computed before starting looping through the triangles */
    vector<Vec3> vert0;
    vector<Vec3> edge1;
    vector<Vec3> edge2;
    vector<Vec3> pvec;
    vector<bool> is_parallel;

    /** Function to find with Moller-Trumbore algorithm whether the ray and the triangle intersect or not */
    const bool ray_intersects_triangle(const Vec3 &origin, const Vec3 &direction, const int face);

    /** Function to reserve memory for precompute data */
    const void reserve(const int n);
};

} /* namespace femocs */

#endif /* TETGENMESH_H_ */