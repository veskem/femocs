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
Pic<dim>::Pic(fch::Laplace<dim> &laplace_solver) : laplace_solver(laplace_solver){
}

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
template<>
int Pic<3>::injectElectrons(const double* const r, const size_t n, FieldReader &fr) {
    for (size_t i = 0; i < n; i++) {
        r_el.push_back(dealii::Point<3>(r[i*3+0],r[i*3+1],r[i*3+2]));
        v_el.push_back(dealii::Point<3>(0.0,0.0,0.0));
        cout << "locating the cell" << endl;
        cid_el.push_back(fr.update_point_cell(r_el[i], 10));
    }
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
void Pic<dim>::pushParticles(const double dt, FieldReader &fr) {

    for (size_t i = 0; i < r_el.size(); i++) {
        //Leapfrog method:
        // positions defined ON the time steps, velocities defined at half time steps
        dealii::Tensor<1,dim> Efield = laplace_solver.probe_efield(r_el[i], cid_el[i]) ; // Get the field!

        v_el[i] = v_el[i] + q_over_m_factor*Efield*(dt*1e15);
        r_el[i] = r_el[i] + v_el[i]*dt;

        //Update the cid_el && check if any particles have left the domain
        cid_el[i] = fr.update_point_cell(r_el[i], cid_el[i]);
    }
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
    out << "Interpolator properties=id:I:1:pos:R:3" << endl;

    for (int i = 0; i < r_el.size(); ++i)
        out << i << " " << r_el[i][0] << " " << r_el[i][1] << " " << r_el[i][2] << endl;

    out.close();
}


//Tell the compiler which types to actually compile, so that they are available for the linker
//template class Pic<2>;
template class Pic<3>;
}

