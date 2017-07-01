/*
 * currents_and_heating.cc
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>	// UMFpack

#include <cassert>
#include <algorithm>

#include "currents_and_heating.h"
#include "utility.h"

namespace fch {
using namespace dealii;

template<int dim>
CurrentsAndHeating<dim>::CurrentsAndHeating(PhysicalQuantities *pq_) :
        fe_current(currents_degree), dof_handler_current(triangulation),
        fe_heat(heating_degree), dof_handler_heat(triangulation), pq(pq_) {
}

template<int dim>
void CurrentsAndHeating<dim>::import_mesh_from_file(const std::string file_name) {
    MeshPreparer<dim> mesh_preparer;

    mesh_preparer.import_mesh_from_file(&triangulation, file_name);
    mesh_preparer.mark_copper_boundary(&triangulation);
}

template<int dim>
void CurrentsAndHeating<dim>::setup_current_system() {

    dof_handler_current.distribute_dofs(fe_current);
    //std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
    //      << std::endl;

    DynamicSparsityPattern dsp(dof_handler_current.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_current, dsp);
    sparsity_pattern_current.copy_from(dsp);

    system_matrix_current.reinit(sparsity_pattern_current);

    solution_current.reinit(dof_handler_current.n_dofs());
    system_rhs_current.reinit(dof_handler_current.n_dofs());
}

template<int dim>
void CurrentsAndHeating<dim>::setup_heating_system() {

    dof_handler_heat.distribute_dofs(fe_heat);
    //std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
    //      << std::endl;

    DynamicSparsityPattern dsp(dof_handler_heat.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_heat, dsp);
    sparsity_pattern_heat.copy_from(dsp);

    system_matrix_heat.reinit(sparsity_pattern_heat);

    solution_heat.reinit(dof_handler_heat.n_dofs());
    system_rhs_heat.reinit(dof_handler_heat.n_dofs());

    // Initialize the solution to ambient temperature
    // This loop is hopefully optimized by the compiler
    for (std::size_t i = 0; i < solution_heat.size(); i++) {
        solution_heat[i] = ambient_temperature;
    }
}


template<int dim>
void CurrentsAndHeating<dim>::assemble_current_system() {
    QGauss<dim> quadrature_formula(quadrature_degree);
    QGauss<dim-1> face_quadrature_formula(quadrature_degree);

    // Current finite element values
    FEValues<dim> fe_values(fe_current, quadrature_formula,
            update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe_current, face_quadrature_formula,
                update_values | update_quadrature_points | update_JxW_values);

    // Temperature finite element values (only for accessing previous iteration solution)
    FEValues<dim> fe_values_heat(fe_heat, quadrature_formula,
            update_values | update_quadrature_points);
    FEFaceValues<dim> fe_face_values_heat(fe_heat, face_quadrature_formula,
                update_values | update_quadrature_points);

    const unsigned int dofs_per_cell = fe_current.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // ---------------------------------------------------------------------------------------------
    // The previous solution values in the cell quadrature points
    std::vector<double> prev_sol_temperature_values(n_q_points);
    // The previous solution values in the face quadrature points
    std::vector<double> prev_sol_face_temperature_values(n_face_q_points);
    // ---------------------------------------------------------------------------------------------

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_current.begin_active(),
            endc = dof_handler_current.end();

    for (; cell != endc; ++cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values_heat.reinit(cell);
        fe_values_heat.get_function_values(solution_heat, prev_sol_temperature_values);

        // ----------------------------------------------------------------------------------------
        // Local matrix assembly
        // ----------------------------------------------------------------------------------------
        for (unsigned int q = 0; q < n_q_points; ++q) {

            double temperature = prev_sol_temperature_values[q];
            double sigma = pq->sigma(temperature);

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
                            * sigma * fe_values.JxW(q));
            }
        }
        // ----------------------------------------------------------------------------------------
        // Local right-hand side assembly
        // ----------------------------------------------------------------------------------------
        // Emission current BC at the copper surface
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                fe_face_values.reinit(cell, f);

                fe_face_values_heat.reinit(cell, f);
                fe_face_values_heat.get_function_values(solution_heat, prev_sol_face_temperature_values);

                double e_field = 5.0;

                for (unsigned int q = 0; q < n_face_q_points; ++q) {

                    double temperature = prev_sol_face_temperature_values[q];
                    double emission_current = pq->emission_current(e_field, temperature);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += (fe_face_values.shape_value(i, q)
                                * emission_current * fe_face_values.JxW(q));
                    }
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix_current.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

            system_rhs_current(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler_current, BoundaryId::copper_bottom,
            ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, system_matrix_current, solution_current, system_rhs_current);
}


template<int dim>
void CurrentsAndHeating<dim>::solve_current(int max_iter, double tol, bool pc_ssor, double ssor_param) {
    SolverControl solver_control(max_iter, tol);
    SolverCG<> solver(solver_control);

    if (pc_ssor) {
        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(system_matrix_current, ssor_param);
        solver.solve(system_matrix_current, solution_current, system_rhs_current, preconditioner);
    } else {
        solver.solve(system_matrix_current, solution_current, system_rhs_current, PreconditionIdentity());
    }
    //std::cout << "   " << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;
}


template<int dim>
void CurrentsAndHeating<dim>::output_results_current(const std::string filename) const {

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler_current);
    data_out.add_data_vector(solution_current, "potential");

    data_out.build_patches();

    try {
        std::ofstream output(filename);
        data_out.write_vtk(output);
    } catch (...) {
        std::cerr << "WARNING: Couldn't open " + filename << ". ";
        std::cerr << "Output is not saved." << std::endl;
    }
}


template class CurrentsAndHeating<2> ;
template class CurrentsAndHeating<3> ;

} // namespace fch

