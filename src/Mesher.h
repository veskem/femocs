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

namespace femocs {
class Bulk;
class Mesh;
class Vacuum;
} /* namespace femocs */

using namespace std;

namespace femocs {

class Mesher {
public:
    Mesher(string mesher, const double latconst);
    virtual ~Mesher() {
    }
    ;

    const void get_test_mesh(Mesh* new_mesh);
    const void get_volume_mesh(Mesh* new_mesh, Bulk* bulk, Vacuum* vacuum, const string cmd);
    const void generate_surf_faces(Mesh* mesh, const int nmax);

    const void separate_vacuum_mesh(Mesh* vacuum_mesh, Mesh* big_mesh, const int nmax,
            const int n_surf, const double zmin, const string cmd);
    const void separate_bulk_mesh(Mesh* bulk_mesh, Mesh* big_mesh, const int nmax, const int n_surf,
            const double zmin, const string cmd);
    const void separate_meshes(Mesh* bulk_mesh, Mesh* vacuum_mesh, Mesh* big_mesh, const int nmax,
            const int n_surf, const double zmin, const string cmd);

    void clean_faces(Mesh* mesh, const double rmax, const string cmd);
    void clean_elems(Mesh* mesh, const double rmax, const string cmd);

    void mark_faces(Mesh* mesh, const AtomReader::Sizes* sizes, const AtomReader::Types* types);
    void mark_faces_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types);
    void mark_elems_bynode(Mesh* mesh, const int nmax, const AtomReader::Types* types);

private:
    const int n_nodes_per_elem = 4;
    const int n_nodes_per_face = 3;
    double latconst;

    const void extract_mesh(vector<bool>* is_vacuum, Mesh* big_mesh, const int nmax,
            const int n_surf, const double zmin);
    void update_list(int* new_list, int* old_list, vector<bool> is_quality, int M);

};

} /* namespace femocs */

#endif /* MESHER_H_ */
