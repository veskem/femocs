/*
 * ProjectHeat.h
 *
 *  Created on: 4.5.2018
 *      Author: Andreas Kyritsakis
 */

#ifndef SRC_PROJECTSPACECHARGE_H_
#define SRC_PROJECTSPACECHARGE_H_

#include "ProjectRunaway.h"
#include "Surface.h"
#include "Interpolator.h"
#include "PoissonSolver.h"
#include "Pic.h"
#include "omp.h"


using namespace std;
namespace femocs {

class ProjectSpaceCharge : public ProjectRunaway {
public:
    ProjectSpaceCharge(AtomReader &reader, Config &conf);
    ~ProjectSpaceCharge(){};

    int run(int timestep);

private:

    /** Run the pic Space charge until convergence is reached */
    int converge_pic();

    /** Find the omega_SC that minimizes the error of factors-currents curve */
    double find_Veff(vector<double> currents);

     void get_currents(double Vappl, vector<double> &curs);

     double get_current_error(vector<double> I_calc, vector<double> I_target);

};

} /* namespace femocs */

#endif /* SRC_PROJECTSPACECHARGE_H_ */
