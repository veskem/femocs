/*
 * Mesher.h
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#ifndef MESHER_H_
#define MESHER_H_

#include <memory>
#include <string>
#include <vector>

#include "../lib/tetgen.h"
#include "Femocs.h"
#include "Surface.h"
#include "Mesh.h"

namespace femocs {
class Vacuum;
} /* namespace femocs */

using namespace std;

namespace femocs {

class Mesher {
public:
    Mesher(string mesher);
    virtual ~Mesher() {
    }
    ;
    const shared_ptr<Mesh> get_union_mesh(shared_ptr<Mesh> mesh_bulk, shared_ptr<Mesh> mesh_volume, const Femocs::SimuCell* cell);
    const shared_ptr<Mesh> get_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum, const string cmd);
    const shared_ptr<Mesh> get_bulk_mesh(shared_ptr<Surface> bulk, const string cmd);

    void clean_faces(shared_ptr<Mesh> mesh, const double rmax, const string cmd);
    void clean_elems(shared_ptr<Mesh> mesh, const double rmax, const string cmd);
    void clean_facets(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell, const string cmd);

    void mark_faces(shared_ptr<Mesh> mesh, shared_ptr<Surface> surf, const Femocs::SimuCell* cell);
    void mark_elems(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell);
    void mark_elems_byvol(shared_ptr<Mesh> mesh, const Femocs::SimuCell* cell);

    void calc_statistics(shared_ptr<Mesh> mesh);

private:
    bool on_face(const double f1n1, const double f1n2, const double f1n3, const double f2);
    bool on_surface(const double f1n1, const double f1n2, const double f1n3, shared_ptr<Surface> surf, const int xn);
    void update_list(int* new_list, int* old_list, vector<bool> is_quality, int M);

};

} /* namespace femocs */

#endif /* MESHER_H_ */
