/*
 * GeneralProject.cpp
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#include "GeneralProject.h"

namespace femocs {

GeneralProject::GeneralProject(const AtomReader &r, const Config& c) :
        reader(r), conf(c), new_mesh(&mesh1), mesh(NULL) {}

} /* namespace femocs */
