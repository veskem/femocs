/*
 * GeneralProject.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include "GeneralProject.h"

namespace femocs {

GeneralProject::GeneralProject() :
        reader(dummy_reader), conf(dummy_config),
        mesh1(&conf.mesh), mesh2(&conf.mesh), new_mesh(&mesh1), mesh(NULL) {}

GeneralProject::GeneralProject(AtomReader &r, Config& c) :
        reader(r), conf(c), mesh1(&c.mesh), mesh2(&c.mesh), new_mesh(&mesh1), mesh(NULL) {}

int GeneralProject::export_data(const double** data, const string& data_type) const {
    if (data_type == LABELS.nodes) {
        *data = mesh->nodes.get();
        return mesh->nodes.size();
    }

    return -1;
}

int GeneralProject::export_data(const int** data, const string& data_type) const {
    if (data_type == LABELS.edges) {
        *data = mesh->edges.get();
        return mesh->edges.size();
    }

    if (data_type == LABELS.triangles) {
        *data = mesh->tris.get();
        return mesh->tris.size();
    }

    if (data_type == LABELS.tetrahedra) {
        *data = mesh->tets.get();
        return mesh->tets.size();
    }

    if (data_type == LABELS.quadrangles) {
        *data = mesh->quads.get();
        return mesh->quads.size();
    }

    if (data_type == LABELS.hexahedra) {
        *data = mesh->hexs.get();
        return mesh->hexs.size();
    }

    return -1;
}

} /* namespace femocs */
