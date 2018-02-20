/*
 * currents_and_heating.cc -> CurrentsAndHeating.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

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

#include <cassert>
#include <algorithm>
#include <cstdlib>

#include "CurrentsAndHeating.h"
#include "Utility.h"

namespace fch {
using namespace dealii;
using namespace std;

// ----------------------------------------------------------------------------------------
/* Class for outputting the current density distribution
 * being calculated from current potential distribution */
template <int dim>
class CurrentPostProcessor : public DataPostprocessorVector<dim> {
public:
    CurrentPostProcessor() : DataPostprocessorVector<dim>("current_density", update_gradients) {}

    void
    compute_derived_quantities_scalar (
            const vector<double>               &/*uh*/,
            const vector<Tensor<1,dim> >       &duh,
            const vector<Tensor<2,dim> >       &/*dduh*/,
            const vector<Point<dim> >          &/*normals*/,
            const vector<Point<dim> >          &/*evaluation_points*/,
            vector<Vector<double> >            &computed_quantities) const {
        for (unsigned int i=0; i<computed_quantities.size(); i++) {
            for (unsigned int d=0; d<dim; ++d)
                computed_quantities[i](d) = duh[i][d];
        }
    }
};
// ----------------------------------------------------------------------------------------
// Class for outputting the electrical conductivity distribution
template <int dim>
class SigmaPostProcessor : public DataPostprocessorScalar<dim> {
    PhysicalQuantities *pq;
public:
    SigmaPostProcessor(PhysicalQuantities *pq_) :
            DataPostprocessorScalar<dim>("conductivity", update_values), pq(pq_) {}

    void
    compute_derived_quantities_scalar (
            const vector<double>               &uh,
            const vector<Tensor<1,dim> >       &/*duh*/,
            const vector<Tensor<2,dim> >       &/*dduh*/,
            const vector<Point<dim> >          &/*normals*/,
            const vector<Point<dim> >          &/*evaluation_points*/,
            vector<Vector<double> >            &computed_quantities) const {
        for (unsigned int i=0; i<computed_quantities.size(); i++) {
            double temperature = uh[i];
            computed_quantities[i](0) = pq->sigma(temperature);
        }
    }
};
// ----------------------------------------------------------------------------

/* ==================================================================
 *  ======================== EmissionSolver ========================
 * ================================================================== */

template<int dim>
EmissionSolver<dim>::EmissionSolver() :
        DealSolver<dim>(), pq(NULL), conf(NULL)
        {}

template<int dim>
EmissionSolver<dim>::EmissionSolver(Triangulation<dim> *tria) :
        DealSolver<dim>(tria), pq(NULL), conf(NULL)
        {}

template<int dim>
void EmissionSolver<dim>::set_bc(const vector<double> &emission) {
    require(emission.size() > 0, "Boundary condition vector should not be empty!");
    bc_values = emission;
    return;
}

template<int dim>
void EmissionSolver<dim>::set_dependencies(PhysicalQuantities *pq_, const femocs::Config::Heating *conf_) {
    pq = pq_;
    conf = conf_;
}

/* ==================================================================
 *  ========================== HeatSolver ==========================
 * ================================================================== */

template<int dim>
HeatSolver<dim>::HeatSolver() :
        EmissionSolver<dim>(), current_solver(NULL) {}

template<int dim>
HeatSolver<dim>::HeatSolver(Triangulation<dim> *tria, const CurrentSolver<dim> *cs) :
        EmissionSolver<dim>(tria), current_solver(cs) {}

template<int dim>
void HeatSolver<dim>::write_vtk(ofstream& out) const {
    SigmaPostProcessor<dim> post_processor(this->pq);
    DataOut<dim> data_out;

    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "temperature");
    data_out.add_data_vector(this->solution, post_processor);

    data_out.build_patches();
    data_out.write_vtk(out);
}

template<int dim>
void HeatSolver<dim>::assemble(const double delta_time) {
    if (this->conf->assemble_method == "euler") {
        assemble_euler_implicit(delta_time);
        this->assemble_rhs(BoundaryId::copper_surface);
    } else
        assemble_crank_nicolson(delta_time);

    this->append_dirichlet(BoundaryId::copper_bottom, this->dirichlet_bc_value);
    this->apply_dirichlet();
}

template<int dim>
void HeatSolver<dim>::assemble_crank_nicolson(const double delta_time) {
    require(false, "Implementation of Crank-Nicolson assembly not verified!");

    require(current_solver, "NULL current solver can't be used!");

    const double gamma = cu_rho_cp / delta_time;
    this->system_matrix = 0;
    this->system_rhs = 0;

    QGauss<dim> quadrature_formula(this->quadrature_degree);
    QGauss<dim-1> face_quadrature_formula(this->quadrature_degree);

    // Heating finite element values
    FEValues<dim> fe_values(this->fe, quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(this->fe, face_quadrature_formula,
            update_values | update_quadrature_points | update_JxW_values);

    // Finite element values for accessing current calculation
    FEValues<dim> fe_values_current(current_solver->fe, quadrature_formula, update_gradients);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // ---------------------------------------------------------------------------------------------
    // The other solution values in the cell quadrature points
    vector<Tensor<1, dim>> potential_gradients(n_q_points);
    vector<Tensor<1, dim>> prev_sol_potential_gradients(n_q_points);
    vector<double> prev_sol_temperature_values(n_q_points);
    vector<Tensor<1, dim>> prev_sol_temperature_gradients(n_q_points);
    // ---------------------------------------------------------------------------------------------

    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator current_cell = current_solver->dof_handler.begin_active();

    int face_index = 0;
    for (; cell != this->dof_handler.end(); ++cell, ++current_cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.get_function_values(this->solution_save, prev_sol_temperature_values);
        fe_values.get_function_gradients(this->solution_save, prev_sol_temperature_gradients);

        fe_values_current.reinit(current_cell);
        fe_values_current.get_function_gradients(current_solver->solution, potential_gradients);
        fe_values_current.get_function_gradients(current_solver->solution_save, prev_sol_potential_gradients);

        // ----------------------------------------------------------------------------------------
        // Local matrix assembly
        // ----------------------------------------------------------------------------------------
        for (unsigned int q = 0; q < n_q_points; ++q) {

            double prev_temperature = prev_sol_temperature_values[q];
            double kappa = this->pq->kappa(prev_temperature);
            double sigma = this->pq->sigma(prev_temperature);

            Tensor<1, dim> prev_temperature_grad = prev_sol_temperature_gradients[q];

            double pot_grad_squared = potential_gradients[q].norm_square();
            double prev_pot_grad_squared = prev_sol_potential_gradients[q].norm_square();

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    cell_matrix(i, j) += fe_values.JxW(q) * (
                            2*gamma*fe_values.shape_value(i, q) * fe_values.shape_value(j, q) // Mass matrix
                            + kappa*fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q) );
                }
                cell_rhs(i) += fe_values.JxW(q) * (
                        2*gamma*fe_values.shape_value(i, q)*prev_temperature
                        - kappa*fe_values.shape_grad(i, q)*prev_temperature_grad
                        + fe_values.shape_value(i, q)*sigma*(pot_grad_squared+prev_pot_grad_squared) );
            }
        }
        // ----------------------------------------------------------------------------------------
        // Local right-hand side assembly
        // ----------------------------------------------------------------------------------------
        // Nottingham BC at the copper surface
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                fe_face_values.reinit(cell, f);
                double nottingham_heat = this->bc_values[face_index++];

                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += (fe_face_values.shape_value(i, q)
                                * 2.0 * nottingham_heat * fe_face_values.JxW(q));
                    }
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                this->system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

            this->system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
}

template<int dim>
void HeatSolver<dim>::assemble_euler_implicit(const double delta_time) {
    require(current_solver, "NULL current solver can't be used!");

    const double gamma = cu_rho_cp / delta_time;
    this->system_matrix = 0;
    this->system_rhs = 0;

    QGauss<dim> quadrature_formula(this->quadrature_degree);

    // Heating finite element values
    FEValues<dim> fe_values(this->fe, quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Finite element values for accessing current calculation
    FEValues<dim> fe_values_current(current_solver->fe, quadrature_formula, update_gradients);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // The other solution values in the cell quadrature points
    vector<Tensor<1, dim>> potential_gradients(n_q_points);
    vector<double> prev_temperatures(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator current_cell = current_solver->dof_handler.begin_active();

    for (; cell != this->dof_handler.end(); ++cell, ++current_cell) {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->solution_save, prev_temperatures);

        // Local matrix assembly
        cell_matrix = 0;
        for (unsigned int q = 0; q < n_q_points; ++q) {
            double temperature = prev_temperatures[q];
            double kappa = this->pq->kappa(temperature);

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    cell_matrix(i, j) += fe_values.JxW(q) * (
                            gamma * fe_values.shape_value(i, q) * fe_values.shape_value(j, q) // Mass matrix
                            + kappa * fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q) );
                }
            }
        }

        fe_values_current.reinit(current_cell);
        fe_values_current.get_function_gradients(current_solver->solution, potential_gradients);

        // Local right-hand-side vector assembly
        cell_rhs = 0;
        for (unsigned int q = 0; q < n_q_points; ++q) {
            double pot_grad_squared = potential_gradients[q].norm_square();
            double temperature = prev_temperatures[q];
            double sigma = this->pq->sigma(temperature);

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                cell_rhs(i) += fe_values.JxW(q) * fe_values.shape_value(i, q)
                        * (gamma * temperature + sigma * pot_grad_squared);
            }
        }

        // Update global matrix and right-hand-side vector
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            this->system_rhs(local_dof_indices[i]) += cell_rhs(i);

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                this->system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
}

template<int dim>
double HeatSolver<dim>::get_face_bc(const unsigned int face) const {
    require(face < this->bc_values.size(), "Invalid index: " + to_string(face));
    return this->bc_values[face];
}

/* ==================================================================
 *  ========================= CurrentSolver ========================
 * ================================================================== */

template<int dim>
CurrentSolver<dim>::CurrentSolver() :
        EmissionSolver<dim>(), heat_solver(NULL) {}

template<int dim>
CurrentSolver<dim>::CurrentSolver(Triangulation<dim> *tria, const HeatSolver<dim> *hs) :
        EmissionSolver<dim>(tria), heat_solver(hs) {}

template<int dim>
void CurrentSolver<dim>::write_vtk(ofstream& out) const {
    CurrentPostProcessor<dim> post_processor; // needs to be before data_out
    DataOut<dim> data_out;

    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "current_potential");
    data_out.add_data_vector(this->solution, post_processor);

    data_out.build_patches();
    data_out.write_vtk(out);
}

template<int dim>
void CurrentSolver<dim>::assemble() {
    assemble_lhs();
    this->assemble_rhs(BoundaryId::copper_surface);
    this->append_dirichlet(BoundaryId::copper_bottom, this->dirichlet_bc_value);
    this->apply_dirichlet();
}

template<int dim>
void CurrentSolver<dim>::assemble_lhs() {
    require(heat_solver, "NULL heat solver can't be used!");

    this->system_matrix = 0;
    this->system_rhs = 0;

    QGauss<dim> quadrature_formula(this->quadrature_degree);

    // Current finite element values
    FEValues<dim> fe_values(this->fe, quadrature_formula,
            update_gradients | update_quadrature_points | update_JxW_values);

    // Temperature finite element values (only for accessing previous iteration solution)
    FEValues<dim> fe_values_heat(heat_solver->fe, quadrature_formula, update_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // The previous temperature values in the cell quadrature points
    vector<double> prev_temperatures(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator heat_cell = heat_solver->dof_handler.begin_active();

    int face_index = 0;
    for (; cell != this->dof_handler.end(); ++cell, ++heat_cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;

        fe_values_heat.reinit(heat_cell);
        fe_values_heat.get_function_values(heat_solver->solution, prev_temperatures);

        // ----------------------------------------------------------------------------------------
        // Local matrix assembly
        // ----------------------------------------------------------------------------------------
        for (unsigned int q = 0; q < n_q_points; ++q) {
            double temperature = prev_temperatures[q];
            double sigma = this->pq->sigma(temperature);

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
                            * sigma * fe_values.JxW(q);
            }
        }
        
        // ----------------------------------------------------------------------------------------
        // Global matrix update
        // ----------------------------------------------------------------------------------------
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                this->system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
}

template<int dim>
double CurrentSolver<dim>::get_face_bc(const unsigned int face) const {
    require(face < this->bc_values.size(), "Invalid index: " + to_string(face));
    return this->bc_values[face];
}

/* ==================================================================
 *  ======================= CurrentHeatSolver ======================
 * ================================================================== */

template<int dim>
CurrentHeatSolver<dim>::CurrentHeatSolver() :
        DealSolver<dim>(), pq(NULL), conf(NULL),
        heat(&this->triangulation, &current),
        current(&this->triangulation, &heat)
{}

template<int dim>
CurrentHeatSolver<dim>::CurrentHeatSolver(PhysicalQuantities *pq_, const femocs::Config::Heating *conf_) :
        DealSolver<dim>(), pq(pq_), conf(conf_),
        heat(&this->triangulation, &current),
        current(&this->triangulation, &heat)
{
    heat.set_dependencies(pq_, conf_);
    current.set_dependencies(pq_, conf_);
}

template<int dim>
void CurrentHeatSolver<dim>::setup(const double temperature) {
    heat.dirichlet_bc_value = temperature;
    current.setup_system();
    heat.setup_system();
}

template<int dim>
void CurrentHeatSolver<dim>::set_dependencies(PhysicalQuantities *pq_, const femocs::Config::Heating *conf_) {
    pq = pq_;
    conf = conf_;
    heat.set_dependencies(pq_, conf_);
    current.set_dependencies(pq_, conf_);
}

template<int dim>
vector<double> CurrentHeatSolver<dim>::get_temperature(const vector<int> &cell_indexes,
        const vector<int> &vert_indexes)
{
    // Initialize vector with a value that is immediately visible if it's not changed to proper one
    vector<double> temperatures(cell_indexes.size(), 1e15);

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&this->triangulation, 0,
                cell_indexes[i], &heat.dof_handler);

        double temperature = heat.solution[dof_cell->vertex_dof_index(vert_indexes[i], 0)];
        temperatures[i] = temperature;
    }
    return temperatures;
}

template<int dim>
vector<Tensor<1, dim>> CurrentHeatSolver<dim>::get_current(const vector<int> &cell_indexes,
        const vector<int> &vert_indexes)
{
    QGauss<dim> quadrature_formula(this->quadrature_degree);
    FEValues<dim> fe_values(current.fe, quadrature_formula, update_gradients);

    vector<Tensor<1, dim> > potential_gradients(quadrature_formula.size());
    const FEValuesExtractors::Scalar potential(0);

    vector<Tensor<1, dim> > currents(cell_indexes.size());

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&this->triangulation, 0,
                cell_indexes[i], &current.dof_handler);

        fe_values.reinit(dof_cell);
        fe_values.get_function_gradients(current.solution, potential_gradients);

        double temperature = heat.solution[dof_cell->vertex_dof_index(vert_indexes[i], 0)];
        Tensor<1, dim> field = -1.0 * potential_gradients.at(vert_indexes[i]);
        Tensor<1, dim> current = pq->sigma(temperature) * field;

        currents[i] = current;
    }
    return currents;
}

template<int dim>
void CurrentHeatSolver<dim>::mark_boundary() {
    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.mark_copper_boundary(&this->triangulation);
}

/* ==================================================================
 * Declare above classes with desired dimensions
 * ================================================================== */

template class EmissionSolver<3> ;
template class HeatSolver<3> ;
template class CurrentSolver<3> ;
template class CurrentHeatSolver<3> ;

} // namespace fch

