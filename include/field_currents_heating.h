/*
 * field_currents_heating.h
 *
 *  Created on: Aug 3, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_FIELD_CURRENTS_HEATING_H_
#define INCLUDE_FIELD_CURRENTS_HEATING_H_


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "physical_quantities.h"

namespace fch {
	using namespace dealii;

	template <int dim>
	class FieldCurrentsHeating {
	public:
		FieldCurrentsHeating (PhysicalQuantities pq_);
		void run();

		Triangulation<dim>* getp_triangulation();

		enum BoundaryId {
			vacuum_top = 2,
			bulk_bottom = 3
		};

		enum MaterialId {
			vacuum_domain = 10,
			copper_domain = 20
		};

	private:

		static bool cell_is_in_vacuum_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell);
		static bool cell_is_in_copper_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell);

		void set_active_fe_indices ();
		void setup_dofs ();
		void assemble_system_newton (const bool first_iteration);
		void solve ();
		void output_results (const unsigned int iteration) const;

		static const unsigned int    field_degree		= 1;
		static const unsigned int    currents_degree	= 1;
		static const unsigned int    heating_degree 	= 1;

		static constexpr double applied_field = 4.0;

		Triangulation<dim>    triangulation;
		FESystem<dim>         vacuum_fe;
		FESystem<dim>         copper_fe;
		hp::FECollection<dim> fe_collection;
		hp::DoFHandler<dim>   dof_handler;

		ConstraintMatrix      constraints;

		SparsityPattern       sparsity_pattern;
		SparseMatrix<double>  system_matrix;

		Vector<double>        solution;
		Vector<double> 		  previous_solution;
		Vector<double>        system_rhs;

		PhysicalQuantities pq;
	};

} // namespace fch

#endif /* INCLUDE_FIELD_CURRENTS_HEATING_H_ */
