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
    vector<double> z_line, V_line, rho_line;

    Point3 apex, anode;
    bool apex_found = false;

    double Vbase, Ebase; //< Values of the input script applied voltage and field before factor multiplication
    double Fmax_base; //< Apex max field for Laplace and apply_factor = 1

    /** Gradually ramp up the field to avoid discontinuities */
    double ramp_field(double start_factor, double target_factor);

    int prepare(int i);

    /** Run the pic Space charge until convergence is reached */
    int converge_pic();

    /** Find the omega_SC that minimizes the error of factors-currents curve */
    double find_Veff();

    void get_currents(double Vappl);

    void prepare_emission();

    double get_current_error();

    void write_output(double Veff);

    void write_line(string filename);

    void full_emission_curve(double Veff);

    void export_on_line();

    void find_apex();

    void write_emission_stats(string filename, bool first_time, double factor);

    void write_emission_data(string filename, bool first_time);

};

} /* namespace femocs */

#endif /* SRC_PROJECTSPACECHARGE_H_ */
