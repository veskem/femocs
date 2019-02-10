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

} /* namespace femocs */
