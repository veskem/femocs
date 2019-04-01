/*
 * ProjectHeat.h
 *
 *  Created on: 4.5.2018
 *      Author: kyritsak
 */

#ifndef SRC_PROJECTHEAT_H_
#define SRC_PROJECTHEAT_H_

#include "ProjectSpaceCharge.h"
#include "Surface.h"
#include "Interpolator.h"
#include "PhysicalQuantities.h"
#include "PoissonSolver.h"
#include "CurrentHeatSolver.h"
#include "Pic.h"
#include "omp.h"

using namespace std;
namespace femocs {

/*
 *
 */
class ProjectHeat : public ProjectSpaceCharge {
public:
    ProjectHeat(AtomReader &reader, Config &conf);
    ~ProjectHeat(){};

    int run(int timestep, double time);

private:

    /** Pick a field solver and calculcate field distribution */
    int run_field_solver();

    /** Pick a heat solver and calculcate temperature & current density distribution */
    int run_heat_solver();

    /** Using constant mesh, solve transient heat and continuity equation until convergence is reached */
    int converge_heat(double T_ambient);

};

} /* namespace femocs */

#endif /* SRC_PROJECTHEAT_H_ */
