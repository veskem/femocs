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

class Mesher {
public:
    Mesher(const string mesher);
    virtual ~Mesher() {
    }
    ;

    const void get_test_mesh(Mesh* new_mesh);
    const void get_volume_mesh(Mesh* new_mesh, Bulk* bulk, Surface* surf, Vacuum* vacuum, const string cmd);
    const void generate_monolayer_surf_faces(Mesh* mesh);
    const void generate_surf_faces(Mesh* mesh);

    const void separate_meshes(Mesh* bulk_mesh, Mesh* vacuum_mesh, Mesh* big_mesh, const int n_bulk,
            const int n_surf, const double zmin, const string cmd);
    const void separate_meshes_bymarker(Mesh* bulk, Mesh* vacuum, Mesh* big_mesh,
            const AtomReader::Types* types, const string cmd);

    const void mark_faces(Mesh* mesh, const AtomReader::Sizes* sizes, const AtomReader::Types* types);
    const void mark_faces_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types);
    const void mark_nodes(Mesh* mesh, const AtomReader::Types* types);
    const void mark_nodes_long(Mesh* mesh, const AtomReader::Types* types);

    const void clean_faces(Mesh* mesh, const double rmax, const string cmd);
    const void clean_elems(Mesh* mesh, const double rmax, const string cmd);

    // Data for precomputation of triangles
    struct PreComp {
        vector<Vec3d> vert0;
        vector<Vec3d> edge1;
        vector<Vec3d> edge2;
        vector<Vec3d> pvec;
        vector<bool> is_parallel;
    };

    PreComp precomp;

private:
    const int n_nodes_per_elem = 4;
    const int n_nodes_per_face = 3;
    const int n_coordinates = 3;
    const double epsilon = 1e-8;

    string mesher;

    const void mark_elems(Mesh* mesh, const AtomReader::Types* types);
    const void post_process_node_marking(Mesh* mesh, const AtomReader::Types* types);

    const vector<bool> get_vacuum_indices(Mesh* big_mesh, const int n_bulk, const int n_surf, const double zmin);

    const void update_list(int* new_list, const int* old_list, const vector<bool> is_quality, const int M);

    const void precompute_triangles(Mesh* mesh, const Vec3d &direction);
    const bool ray_intersects_triangle(const Vec3d &origin, const Vec3d &direction, const int face);
    const int ray_intersects_surface(Mesh* mesh, const Vec3d &origin, const Vec3d &direction);
    const bool ray_intersects_surface_fast(Mesh* mesh, const Vec3d &origin, const Vec3d &direction);
};

} /* namespace femocs */

#endif /* MESHER_H_ */
