/*
 * Pic.cpp
 *
 *  Created on: 17.01.2018
 *      Author: Kyrre, Andreas
 */

#include "Pic.h"
#include "Macros.h"

#include <deal.II/base/tensor.h>

namespace femocs {
template<int dim>
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver, FieldReader &fr) :
laplace_solver(laplace_solver), fr(fr){}

template<int dim>
Pic<dim>::~Pic() {

}

//template<>
//int Pic<2>::injectElectrons(const double* const r, const size_t n, FieldReader &fr) {
//    for (size_t i = 0; i < n; i++) {
//        r_el.push_back(dealii::Point<2>(r[i*2+0],r[i*2+1]));
//        v_el.push_back(dealii::Point<2>(0.0,0.0));
//        cid_el.push_back(fr.update_point_cell(r_el[i], 10));
//    }
//}
template<int dim>
int Pic<dim>::injectElectrons(const double* const r, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        r_el.push_back(dealii::Point<3>(r[i*3+0],r[i*3+1],r[i*3+2]));
        v_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
        cout << "locating the cell" << endl;
        cid_el.push_back(fr.update_point_cell(r_el[i], 10));
    }
}

template<int dim>
int Pic<dim>::injectElectrons(const fch::CurrentsAndHeating<3> &ch_solver, const double &dt_pic) {

    std::vector<dealii::Point<dim>> new_el = ch_solver.inject_electrons(dt_pic);

    for (auto& point : new_el){
        r_el.push_back(point);
        v_el.push_back(0. * point);

        cout << "locating point ";
        cid_el.push_back(fr.update_point_cell(point, 10));
        cout << cid_el.back() << endl;
        //std::printf("inserting electron in point %e, %e, %e\n", point[0], point[1], point[2]);
    }
    return 0;

}

template<int dim>
void Pic<dim>::computeField(const double E0) {
    //Call the laplace solver with the list of positions and charge(s)
    double t0;

    start_msg(t0, "=== Initializing Laplace solver...");
    laplace_solver.set_applied_efield(-E0);
    laplace_solver.setup_system();
    end_msg(t0);

    start_msg(t0, "Assembling system lhs");
    laplace_solver.assemble_system_lhs();
    end_msg(t0);

    start_msg(t0, "Assembling Neumann boundary rhs");
    laplace_solver.assemble_system_neuman(fch::BoundaryId::vacuum_top);
    end_msg(t0);

    start_msg(t0, "Assembling charge rhs");
    laplace_solver.assemble_system_pointcharge(r_el, -q_over_eps0*Wsp, cid_el);
    end_msg(t0);

    start_msg(t0, "Applying Dirichlet boundary conditions");
    laplace_solver.assemble_system_dirichlet(fch::BoundaryId::copper_surface, 0.0);
    end_msg(t0);

    start_msg(t0, "Closing the FEM system of equations");
    laplace_solver.assemble_system_finalize();
    end_msg(t0);

    start_msg(t0, "Solving Poisson equation");
    laplace_solver.solve();
    end_msg(t0);




}

template<int dim>
void Pic<dim>::pushParticles(const double dt) {

    for (size_t i = 0; i < r_el.size(); i++) {
        //Skip lost particles
        for (auto lost : lost_el) {
            if (lost == i) continue;
        }

        //Leapfrog method:
        // positions defined ON the time steps, velocities defined at half time steps
        cout << "probig field at particle location...\n";
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(r_el[i], cid_el[i]) ; // Get the field!

        v_el[i] = v_el[i] + q_over_m_factor*Efield*(dt);
        r_el[i] = r_el[i] + v_el[i]*dt;

        //Update the cid_el && check if any particles have left the domain
        cout << "updating cell for particles\n";
        cid_el[i] = fr.update_point_cell(r_el[i], cid_el[i]);
        if (cid_el[i] == -1) {
            lost_el.push_back(i);
        }
    }
}

template<int dim>
void Pic<dim>::clearLostParticles(){
    size_t npart = r_el.size();
    size_t nlost = 0;

    //Delete the lost particles from the arrays
    for (size_t i = 0; i < npart; i++) {
        bool islost=false;
        //Is this particle lost?
        for (auto lost : lost_el) {
            if (lost == i) {
                islost=true;
                nlost++;
                break;
            }
        }
        if (nlost==0 or islost) continue; // Don't shuffle this particle left

        r_el[i-nlost] = r_el[i];
        v_el[i-nlost] = v_el[i];
        cid_el[i-nlost] = cid_el[i];
    }

    //Shrink the arrays
    if (nlost > 0){
        r_el.resize(npart-nlost);
        v_el.resize(npart-nlost);
        cid_el.resize(npart-nlost);
        cout << "Particles where lost! nlost=" << nlost << endl;
    }

    lost_el.clear();
}

template<int dim>
void Pic<dim>::writeParticles(const string filename) {
    ofstream out;
    out.setf(std::ios::scientific);
    out.precision(6);

    string ftype = get_file_type(filename);
    if (ftype == "movie") out.open(filename, ios_base::app);
    else out.open(filename);
    require(out.is_open(), "Can't open a file " + filename);

    out << r_el.size() << endl;
    out << "Interpolator properties=id:I:1:pos:R:3:vel:R:3:cell:I:1" << endl;

    for (int i = 0; i < r_el.size(); ++i)
        out << i << " " << r_el[i][0] << " " << r_el[i][1] << " " << r_el[i][2] << " " <<
        v_el[i][0] << " " << v_el[i][1] << " " << v_el[i][2] << " " << cid_el[i] << endl;

    out.close();
}


//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;
}
