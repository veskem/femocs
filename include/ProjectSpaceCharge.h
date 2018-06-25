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

    int run(int timestep, double time);

private:
    vector<double> I_pic, I_sc;

    /** Run the pic Space charge until convergence is reached */
    int converge_pic();

    /** Find the omega_SC that minimizes the error of factors-currents curve */
    double find_Veff();

     void get_currents(double Vappl);

     double get_current_error();

     void write_results(double Veff);

};

} /* namespace femocs */

#endif /* SRC_PROJECTSPACECHARGE_H_ */
