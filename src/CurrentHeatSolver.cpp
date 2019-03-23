/*
 * currents_and_heating.cc -> CurrentsAndHeating.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan, Mihkel
 */

#include <deal.II/numerics/data_out.h>
#include <deal.II/base/work_stream.h>

#include "CurrentHeatSolver.h"
#include "EmissionReader.h"


namespace femocs {
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
        DealSolver<dim>(), pq(NULL), conf(NULL), bc_values(NULL)
        {}

template<int dim>
EmissionSolver<dim>::EmissionSolver(Triangulation<dim> *tria, vector<double>* bcs) :
        DealSolver<dim>(tria), pq(NULL), conf(NULL), bc_values(bcs)
        {}

/* ==================================================================
 *  ========================== HeatSolver ==========================
 * ================================================================== */

template<int dim>
HeatSolver<dim>::HeatSolver() :
        EmissionSolver<dim>(), current_solver(NULL), one_over_delta_time(0) {}

template<int dim>
HeatSolver<dim>::HeatSolver(Triangulation<dim> *tria, const CurrentSolver<dim> *cs, vector<double>* bcs) :
        EmissionSolver<dim>(tria, bcs), current_solver(cs), one_over_delta_time(0) {}

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
void HeatSolver<dim>::setup_system() {
    DealSolver<dim>::setup_system();

    const unsigned int n_dofs = this->size();
    joule_heat.reinit(n_dofs);
    total_heat.reinit(n_dofs);
    this->dof_volume.resize(n_dofs);
}

template<int dim>
void HeatSolver<dim>::assemble(const double delta_time) {
    require(current_solver, "NULL current solver can't be used!");
    require(delta_time > 0, "Invalid delta time: " + d2s(delta_time));

    this->one_over_delta_time = 1.0 / delta_time;
    this->system_matrix = 0;
    this->system_rhs = 0;

    LinearSystem system(&this->system_rhs, &this->system_matrix);
    QGauss<dim> quadrature_formula(this->quadrature_degree);

    const unsigned int n_dofs = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    WorkStream::run(this->dof_handler.begin_active(),this->dof_handler.end(),
            std::bind(&HeatSolver<dim>::assemble_local_cell,
                    this,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3),
            std::bind(&HeatSolver<dim>::copy_global_cell,
                    this,
                    std::placeholders::_1,
                    std::ref(system)),
            ScratchData(this->fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values),
            CopyData(n_dofs, n_q_points)
    );

    if (this->write_time()) {
        this->joule_heat = this->system_rhs;
        this->calc_dof_volumes();
    }
    this->assemble_rhs(BoundaryID::copper_surface);
    if (this->write_time()) this->total_heat = this->system_rhs;
    this->append_dirichlet(BoundaryID::copper_bottom, this->dirichlet_bc_value);
    this->apply_dirichlet();
}

template<int dim>
void HeatSolver<dim>::assemble_crank_nicolson(const double delta_time) {
    require(false, "Implementation of Crank-Nicolson assembly not verified!");

    /* TODO:
     * change temperature_grad to temperature
     * add prev_nottingham values
     * interpolate prev_potential and prev_nottingham for current mesh
     *   prev_potential could be read from heat_transfer
     *   for nottingham something else should be done...
     */

    /*
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
    vector<Tensor<1, dim>> potential_grads(n_q_points);
    vector<Tensor<1, dim>> prev_potential_grads(n_q_points);
    vector<double> prev_temperatures(n_q_points);
    vector<Tensor<1, dim>> prev_temperature_grads(n_q_points);
    // ---------------------------------------------------------------------------------------------

    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator current_cell = current_solver->dof_handler.begin_active();

    int face_index = 0;
    for (; cell != this->dof_handler.end(); ++cell, ++current_cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.get_function_values(this->solution_save, prev_temperatures);
        fe_values.get_function_gradients(this->solution_save, prev_temperature_grads);

        fe_values_current.reinit(current_cell);
        fe_values_current.get_function_gradients(current_solver->solution, potential_grads);
        fe_values_current.get_function_gradients(current_solver->solution_save, prev_potential_grads);

        // ----------------------------------------------------------------------------------------
        // Local matrix assembly
        // ----------------------------------------------------------------------------------------
        for (unsigned int q = 0; q < n_q_points; ++q) {

            double temperature = prev_temperatures[q];
            double kappa = this->pq->kappa(temperature);
            double sigma = this->pq->sigma(temperature);

            Tensor<1, dim> temperature_grad = prev_temperature_grads[q];

            double pot_grad_squared = potential_grads[q].norm_square();
            double prev_pot_grad_squared = prev_potential_grads[q].norm_square();

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    cell_matrix(i, j) += fe_values.JxW(q) * (
                            gamma*fe_values.shape_value(i, q) * fe_values.shape_value(j, q) // Mass matrix
                            + 0.5*kappa*fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q) );
                }
                cell_rhs(i) += fe_values.JxW(q) * (
                        gamma*fe_values.shape_value(i, q)*temperature
                        - 0.5*kappa*fe_values.shape_grad(i, q)*temperature_grad  // TODO check this
                        + 0.5*sigma*fe_values.shape_value(i, q)*(pot_grad_squared+prev_pot_grad_squared) );
            }
        }

        // ----------------------------------------------------------------------------------------
        // Local right-hand side assembly
        // ----------------------------------------------------------------------------------------
        // Nottingham BC at the copper surface
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->boundary_id() == BoundaryID::copper_surface) {
                fe_face_values.reinit(cell, f);
                double nottingham_heat = this->get_face_bc(face_index++);

                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += fe_face_values.shape_value(i, q)
                                * nottingham_heat * fe_face_values.JxW(q);
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
    //*/
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

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // The other solution values in the cell quadrature points
    vector<Tensor<1, dim>> potential_gradients(n_q_points);
    vector<double> prev_temperatures(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();

    for (; cell != this->dof_handler.end(); ++cell) {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->solution, prev_temperatures);
        fe_values.get_function_gradients(current_solver->solution, potential_gradients);

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
void HeatSolver<dim>::assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
        ScratchData &scratch_data, CopyData &copy_data) const
{
    const unsigned int n_dofs = copy_data.n_dofs;
    const unsigned int n_q_points = copy_data.n_q_points;

    const double gamma = cu_rho_cp * one_over_delta_time;

    // The other solution values in the cell quadrature points
    vector<Tensor<1, dim>> potential_gradients(n_q_points);
    vector<double> prev_temperatures(n_q_points);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values.get_function_values(this->solution, prev_temperatures);
    scratch_data.fe_values.get_function_gradients(current_solver->solution, potential_gradients);

    // Local matrix assembly
    copy_data.cell_matrix = 0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
        double temperature = prev_temperatures[q];
        double kappa = this->pq->kappa(temperature);

        for (unsigned int i = 0; i < n_dofs; ++i) {
            for (unsigned int j = 0; j < n_dofs; ++j) {
                copy_data.cell_matrix(i, j) += scratch_data.fe_values.JxW(q) * (
                        gamma * scratch_data.fe_values.shape_value(i, q) * scratch_data.fe_values.shape_value(j, q) // Mass matrix
                        + kappa * scratch_data.fe_values.shape_grad(i, q) * scratch_data.fe_values.shape_grad(j, q) );
            }
        }
    }

    // Local right-hand-side vector assembly
    copy_data.cell_rhs = 0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
        double pot_grad_squared = potential_gradients[q].norm_square();
        double temperature = prev_temperatures[q];
        double sigma = this->pq->sigma(temperature);

        for (unsigned int i = 0; i < n_dofs; ++i) {
            copy_data.cell_rhs(i) += scratch_data.fe_values.JxW(q) * scratch_data.fe_values.shape_value(i, q)
                    * (gamma * temperature + sigma * pot_grad_squared);
        }
    }

    // Obtain dof indices for updating global matrix and right-hand-side vector
    cell->get_dof_indices(copy_data.dof_indices);
}

/* ==================================================================
 *  ========================= CurrentSolver ========================
 * ================================================================== */

template<int dim>
CurrentSolver<dim>::CurrentSolver() :
        EmissionSolver<dim>(), heat_solver(NULL) {}

template<int dim>
CurrentSolver<dim>::CurrentSolver(Triangulation<dim> *tria, const HeatSolver<dim> *hs, vector<double> *bcs) :
        EmissionSolver<dim>(tria, bcs), heat_solver(hs) {}

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
    require(heat_solver, "NULL heat solver can't be used!");

    this->system_matrix = 0;
    this->system_rhs = 0;

    LinearSystem system(&this->system_rhs, &this->system_matrix);
    QGauss<dim> quadrature_formula(this->quadrature_degree);

    const unsigned int n_dofs = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    WorkStream::run(this->dof_handler.begin_active(),this->dof_handler.end(),
            std::bind(&CurrentSolver<dim>::assemble_local_cell,
                    this,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3),
            std::bind(&CurrentSolver<dim>::copy_global_cell,
                    this,
                    std::placeholders::_1,
                    std::ref(system)),
            ScratchData(this->fe, quadrature_formula, update_gradients | update_quadrature_points | update_JxW_values),
            CopyData(n_dofs, n_q_points)
    );

    this->assemble_rhs(BoundaryID::copper_surface);
    this->append_dirichlet(BoundaryID::copper_bottom, this->dirichlet_bc_value);
    this->apply_dirichlet();
}

template<int dim>
void CurrentSolver<dim>::assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
        ScratchData &scratch_data, CopyData &copy_data) const
{
    const unsigned int n_dofs = copy_data.n_dofs;
    const unsigned int n_q_points = copy_data.n_q_points;

    // The previous temperature values in the cell quadrature points
    vector<double> prev_temperatures(n_q_points);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values.get_function_values(heat_solver->solution, prev_temperatures);

    // Local matrix assembly
    copy_data.cell_matrix = 0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
        double temperature = prev_temperatures[q];
        double sigma = this->pq->sigma(temperature);

        for (unsigned int i = 0; i < n_dofs; ++i) {
            for (unsigned int j = 0; j < n_dofs; ++j) {
                copy_data.cell_matrix(i, j) += sigma * scratch_data.fe_values.JxW(q) *
                scratch_data.fe_values.shape_grad(i, q) * scratch_data.fe_values.shape_grad(j, q);
            }
        }
    }

    // Nothing to add to local right-hand-side vector assembly
//    copy_data.cell_rhs = 0;

    // Obtain dof indices for updating global matrix and right-hand-side vector
    cell->get_dof_indices(copy_data.dof_indices);
}

/* ==================================================================
 *  ======================= CurrentHeatSolver ======================
 * ================================================================== */

template<int dim>
CurrentHeatSolver<dim>::CurrentHeatSolver() :
        DealSolver<dim>(),
        heat(&this->triangulation, &current, NULL),
        current(&this->triangulation, &heat, NULL),
        pq(NULL), conf(NULL)
{}

template<int dim>
CurrentHeatSolver<dim>::CurrentHeatSolver(PhysicalQuantities *pq_, const Config::Heating *conf_, EmissionReader *emission) :
        DealSolver<dim>(),
        heat(&this->triangulation, &current, emission->get_nottingham()),
        current(&this->triangulation, &heat, emission->get_current_densities()),
        pq(pq_), conf(conf_)
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
void CurrentHeatSolver<dim>::set_dependencies(PhysicalQuantities *pq_, const Config::Heating *conf_) {
    pq = pq_;
    conf = conf_;
    heat.set_dependencies(pq_, conf_);
    current.set_dependencies(pq_, conf_);
}

template<int dim>
void CurrentHeatSolver<dim>::export_temp_rho(vector<double> &temp, vector<Tensor<1,dim>> &rho) const {
    heat.export_solution(temp);        // extract temperatures
    current.export_solution_grad(rho); // extract fields

    // transfer fields to current densities
    for (int i = 0; i < temp.size(); ++i)
        rho[i] = pq->sigma(temp[i]) * rho[i];
}

template<int dim>
void CurrentHeatSolver<dim>::mark_mesh() {
    this->mark_boundary(BoundaryID::copper_surface, BoundaryID::copper_bottom,
            BoundaryID::copper_sides, BoundaryID::copper_surface);
}

template<int dim>
void CurrentHeatSolver<dim>::write_xyz(ofstream& out) const {
    // write the start of xyz header
    FileWriter::write_xyz(out);

    // write Ovito header
    out << "properties=id:I:1:pos:R:3:force:R:3:total_heat:R:1:nottingham_heat:R:1"
            ":joule_heat:R:1:volume:R:1:temperature:R:1:potential:R:1:rho:R:1\n";

    // extract coordinates of dofs
    vector<Point<dim>> support_points;
    heat.export_dofs(support_points);

    // extract temperatures and current densities
    vector<double> temp;
    vector<Tensor<1,dim>> rho;
    export_temp_rho(temp, rho);

    const int n_dofs = support_points.size();
    const int n_verts = heat.tria->n_used_vertices();

    // generate dof index -> vertex index mapping
    vector<int> dof2vertex(n_dofs);
    for (int i = 0; i < n_verts; ++i)
        dof2vertex[heat.vertex2dof[i]] = i;

    // write data
    for (int i = 0; i < n_dofs; ++i) {
        out << dof2vertex[i] << " " << Point3(support_points[i]) << " " << Vec3(rho[dof2vertex[i]])
                << " " << heat.total_heat(i) << " " << heat.total_heat(i) - heat.joule_heat(i)
                << " " << heat.joule_heat(i) << " " << heat.dof_volume[i] << " " << heat.solution(i)
                << " " << current.solution(i) << " " << rho[dof2vertex[i]].norm() << "\n";
    }
}

/* ==================================================================
 * Declare above classes with desired dimensions
 * ================================================================== */

template class EmissionSolver<3> ;
template class HeatSolver<3> ;
template class CurrentSolver<3> ;
template class CurrentHeatSolver<3> ;

} // namespace femocs
