/*
 * Mesher.h
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#ifndef MESHER_H_
#define MESHER_H_

#include <TetgenCells.h>
#include "Macros.h"
#include "AtomReader.h"
#include "Primitives.h"
#include "Media.h"
#include "Mesh.h"


using namespace std;
namespace femocs {

/** Class to mark Mesh nodes with ray-triangle intersection technique,
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle
 */
class RaySurfaceIntersect {
public:
    /** Constructor of RaySurfaceIntersect  */
    RaySurfaceIntersect(Mesh* mesh);
    /** Unimplemented destructor of RaySurfaceIntersect */
    virtual ~RaySurfaceIntersect() {};

    /** Function to find with Moller-Trumbore algorithm how many times the ray and the surface intersect */
    const int ray_intersects_surface(const Vec3 &origin, const Vec3 &direction);

    /** Function to find with Moller-Trumbore algorithm whether the ray and the surface intersect or not */
    const bool ray_intersects_surface_fast(const Vec3 &origin, const Vec3 &direction);
    
    /** Function to precompute the data needed to execute the Moller-Trumbore algorithm */
    const void precompute_triangles(const Vec3 &direction);

private:
    /** Constants to specify the tolerances */
    const double epsilon = 1e-3;
    const double zero = -1.0 * epsilon;
    const double one = 1.0 + epsilon;

    /** Pointer to Mesh with nodes and surface faces */
    Mesh* mesh;

    /** Data computed before starting looping through the triangles */
    vector<Vec3> vert0;
    vector<Vec3> edge1;
    vector<Vec3> edge2;
    vector<Vec3> pvec;
    vector<bool> is_parallel;

    /** Corner of mean plane */
    Point3 corner_node;

    /** Function to find with Moller-Trumbore algorithm whether the ray and the triangle intersect or not */
    const bool ray_intersects_triangle(const Vec3 &origin, const Vec3 &direction, const int face);

    /** Function to reserve memory for precompute data */
    const void reserve(const int n);
};


class Mesher {
public:
    Mesher(Mesh* mesh);
    virtual ~Mesher() {
    }
    ;

    /** Function to generate mesh from surface, bulk and vacuum atoms.
     * Return value indicates the success of it. */
    const bool generate_mesh(Bulk &bulk, Surface &surf, Vacuum &vacuum, const string& cmd);

    const void generate_mesh_appendices();

    const bool mark_mesh(const bool postprocess);

    const void separate_meshes(Mesh* vacuum, Mesh* bulk, const string& cmd);
    const void separate_meshes(Mesh* vacuum, const string& cmd);

    const void separate_meshes_vol2(Mesh* vacuum, Mesh* bulk, const string& cmd);

    const void separate_meshes_noclean(Mesh* vacuum, Mesh* bulk, const string& cmd);
    const void separate_meshes_noclean(Mesh* vacuum, const string& cmd);

    const void swap_sharp_elements(vector<bool> &elem_in_vacuum, vector<bool> &elem_on_perim);

private:
    Mesh* mesh;

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

} /* namespace femocs */

#endif /* MESHER_H_ */
