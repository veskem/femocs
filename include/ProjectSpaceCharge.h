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
    vector<double> fact_full, I_full;
    vector<double> z_line, V_line, rho_line;

    Point3 apex, anode;
    bool apex_found = false;

    /** Run the pic Space charge until convergence is reached */
    int converge_pic(bool ramp = false);

    /** Find the omega_SC that minimizes the error of factors-currents curve */
    double find_Veff();

     void get_currents(double Vappl);

     double get_current_error();

     void write_output(double Veff);

     void write_line(string filename);

     void full_curve(double Veff);

     void find_Wsp();

     void export_on_line();

     void find_apex();


};

} /* namespace femocs */

#endif /* SRC_PROJECTSPACECHARGE_H_ */
