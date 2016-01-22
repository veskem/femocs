/*
 * Vacuum.h
 *
 *  Created on: 21.1.2016
 *      Author: veske
 */

#ifndef VACUUM_H_
#define VACUUM_H_

#include <memory>

#include "Femocs.h"
#include "Medium.h"
#include "Surface.h"

using namespace std;
namespace femocs {

/**
 *
 */
class Vacuum: public Medium {
public:
    Vacuum();
    const void generate(const Femocs::SimuCell* cell, shared_ptr<Surface> surf);

private:
    const void initialise_vars(const int natoms);
    const void addpoint(const double x, const double y, const double z, const int type);
};

} /* namespace femocs */

#endif /* VACUUM_H_ */