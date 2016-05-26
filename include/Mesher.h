/*
 * Mesher.h
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#ifndef MESHER_H_
#define MESHER_H_

#include "Macros.h"
#include "AtomReader.h"
#include "Media.h"
#include "Mesh.h"
#include "Primitives.h"

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
    const double epsilon = 1e-8;
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

    /** Function to generate simple mesh that consists of one tetrahedron */
    const void generate_simple(const string cmd);

    /** Function to generate mesh from surface, bulk and vacuum atoms */
    const void generate_mesh(Bulk &bulk, Surface &surf, Vacuum &vacuum, const string cmd);

    const void generate_monolayer_surf_faces();
    const void generate_surf_faces();

    const void separate_meshes_byseq(Mesh* bulk, Mesh* vacuum, const string cmd);
    const void separate_meshes(Mesh* bulk, Mesh* vacuum, const AtomReader::Types* types, const string cmd);

    const void mark_mesh(const AtomReader::Types* types, const bool postprocess);
    const void mark_nodes_long(const AtomReader::Types* types);

private:
    Mesh* mesh;

    const void mark_faces(const AtomReader::Types* types);
    const void mark_elems(const AtomReader::Types* types);
    const void post_process_marking(const AtomReader::Types* types);

    const vector<bool> get_vacuum_indices();
    const void update_list(int* new_list, const int* old_list, const vector<bool> is_quality, const int M);

};

} /* namespace femocs */

#endif /* MESHER_H_ */
