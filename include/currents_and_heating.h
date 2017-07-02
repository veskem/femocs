/*
 * currents_and_heating.h
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_CURRENTS_AND_HEATING_H_
#define INCLUDE_CURRENTS_AND_HEATING_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>
#include <iostream>
#include <tuple>

#include "mesh_preparer.h"
#include "physical_quantities.h"
#include "laplace.h"

namespace fch {

// forward declaration for Laplace to exist when declaring CurrentsAndHeating
template<int dim> class Laplace;

using namespace dealii;


template<int dim>
class CurrentsAndHeating {
public:
    CurrentsAndHeating(double time_step_, PhysicalQuantities *pq_);

    /**
     * Imports mesh from file and sets the boundary indicators corresponding to copper
     * @param file_name file from the mesh is imported
     */
    void import_mesh_from_file(const std::string file_name);

    /** Sets up degrees of freedom and the sparsity pattern */
    void setup_current_system();
    void setup_heating_system();

    void assemble_current_system();
    void assemble_heating_system();

    /** solves the current/heat matrix equation using conjugate gradients
     * @param max_iter maximum number of iterations allowed
     * @param tol tolerance
     * @param pc_ssor flag to use SSOR preconditioner
     * @param ssor_param parameter to SSOR preconditioner, 1.2 is known to work well with laplace
     */
    unsigned int solve_current(int max_iter = 2000, double tol = 1e-9, bool pc_ssor = true,
            double ssor_param = 1.2);
    unsigned int solve_heat(int max_iter = 2000, double tol = 1e-9, bool pc_ssor = true,
            double ssor_param = 1.2);

    /** outputs the results to a specified file */
    void output_results_current(const std::string filename = "current_solution.vtk") const;
    void output_results_heating(const std::string filename = "heating_solution.vtk") const;

    /** Set the electric field boundary condition on copper-vacuum boundary */
    void set_electric_field_bc(const Laplace<dim> &laplace);
    void set_electric_field_bc(const std::vector<double> &e_fields);

    double get_max_temperature();

private:

    static constexpr unsigned int currents_degree = 1;
    static constexpr unsigned int heating_degree = 1;

    static constexpr double ambient_temperature = 300.0;

    static constexpr double cu_rho_cp = 3.4496e-21; // J/(K*nm^3)

    double time_step;

    Triangulation<dim> triangulation;

    // Current specific variables
    FE_Q<dim> fe_current;
    DoFHandler<dim> dof_handler_current;
    SparsityPattern sparsity_pattern_current;
    SparseMatrix<double> system_matrix_current;
    Vector<double> solution_current;
    Vector<double> old_solution_current;
    Vector<double> system_rhs_current;

    // Heating specific variables
    FE_Q<dim> fe_heat;
    DoFHandler<dim> dof_handler_heat;
    SparsityPattern sparsity_pattern_heat;
    SparseMatrix<double> system_matrix_heat;
    Vector<double> solution_heat;
    Vector<double> old_solution_heat;
    Vector<double> system_rhs_heat;

    Vector<double> const_temperature_solution;

    PhysicalQuantities *pq;

    /** Mapping of copper interface face to vacuum side e field norm
     * (copper_cell_index, copper_cell_face) <-> (electric field norm)
     */
    std::map<std::pair<unsigned, unsigned>, double> interface_map_field;
};

} // end fch namespace

#endif /* INCLUDE_CURRENTS_AND_HEATING_H_ */
