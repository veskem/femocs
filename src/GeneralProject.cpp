/*
 * GeneralProject.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include "GeneralProject.h"

namespace femocs {

GeneralProject::GeneralProject() :
        fields(NULL), temperatures(NULL), forces(NULL),
        reader(dummy_reader), conf(dummy_config), new_mesh(&mesh1), mesh(NULL) {}

GeneralProject::GeneralProject(AtomReader &r, Config& c) :
        fields(NULL), temperatures(NULL), forces(NULL),
        reader(r), conf(c), new_mesh(&mesh1), mesh(NULL) {}

} /* namespace femocs */
