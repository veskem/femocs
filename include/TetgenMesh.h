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
 * Class to create and handle finite element mesh in Tetgen format,
 * http://wias-berlin.de/software/tetgen/1.5/ */
class TetgenMesh {
public:
    TetgenMesh() {}
    ~TetgenMesh() {}

    /** Function to generate simple mesh that consists of one element */
    const bool generate_simple();

    /** Function to generate mesh from surface, bulk and vacuum atoms */
    const bool generate(const Media& bulk, const Media& surf, const Media& vacuum, const string& cmd);

    /** Separate tetrahedra into hexahedra by adding node to the centroid of the
     * tetrahedron edges, nodes and tetrahedron itself */
    const bool generate_appendices();

    /** Mark mesh nodes and elements by their location relative to the surface mesh */
    const bool mark_mesh(const bool postprocess);

    /** Separate generated mesh into bulk and vacuum meshes */
    const bool separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd);

    /** Group hexahedra around central tetrahedral node */
    const void group_hexahedra();

    /** Generate list of nodes that surround the tetrahedral nodes */
    const vector<vector<unsigned int>> get_voronoi_cells() const;

    /** Smoothen hexahedra on the surface */
    const bool smoothen(double radius, double smooth_factor, double r_cut);

    /** Use tetgen built-in function to write elements, faces, edges and nodes into file */
    const bool write_tetgen(const string& file_name);

    /** Copy node and element data from Tetgen input buffer into output one */
    const bool recalc();

    /** Perform Tetgen calculation on input buffer and store it in output one */
    const bool recalc(const string& cmd);

    /** Perform double Tetgen calculation on input buffer and store it in output one */
    const bool recalc(const string& cmd1, const string& cmd2);

    /** Objects holding operations for accessing cell data */
    TetgenNodes nodes = TetgenNodes(&tetIOout, &tetIOin);
    TetgenEdges edges = TetgenEdges(&tetIOout);
    TetgenFaces faces = TetgenFaces(&tetIOout);
    TetgenElements elems = TetgenElements(&tetIOout, &tetIOin);
    Hexahedra hexahedra = Hexahedra(&tetIOout);

    const int n_coordinates = 3;     ///< Number of coordinates
    const int n_edges_per_face = 3;  ///< Number of edges on a triangle
    const int n_edges_per_elem = 6;  ///< Number of edges on a tetrahedron
    const int n_faces_per_elem = 4;  ///< Number of triangles on a tetrahedron

    /** String stream prints the statistics about the mesh */
    friend std::ostream& operator <<(std::ostream &s, const TetgenMesh &t) {
        s << "#hexs=" << t.hexahedra.size()
                << ",\t#tets=" << t.elems.size()
                << ",\t#faces=" << t.faces.size()
                << ",\t#edges=" << t.edges.size()
                << ",\t#nodes=" << t.nodes.size();
        return s;
    }

private:
    tetgenio tetIOin;   ///< Writable mesh data in Tetgen format
    tetgenio tetIOout;  ///< Readable mesh data in Tetgen format

    const void mark_nodes();
    const void mark_edges();
    const void mark_faces();
    const void mark_elems();
    const void mark_elems_byrsi();
    const void mark_elems_vol2();

    const void remark_perimeter_nodes();

    const void remark_elems(const int skip_type);

    const bool post_process_marking();

    const void generate_edges();

    /** Generate surface faces from already existing elements and surface nodes */
    const void generate_surf_faces();
};

/** Class to mark mesh nodes with ray-triangle intersection technique,
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle */
class RaySurfaceIntersect {
public:
    /** Constructor of RaySurfaceIntersect  */
    RaySurfaceIntersect(TetgenMesh* mesh);

    /** Function to find with Moller-Trumbore algorithm whether the ray and the surface intersect or not */
    const bool ray_intersects_surface(const Vec3 &origin, const Vec3 &direction);

    /** Function to precompute the data needed to execute the Moller-Trumbore algorithm */
    const void precompute_triangles(const Vec3 &direction);

private:
    /** Constants to specify the tolerances */
    const double epsilon = 1e-3;
    const double zero = 0.0 - epsilon;
    const double one  = 1.0 + epsilon;

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
