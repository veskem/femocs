/*
 * currents_and_heating.h
 *
 *  Created on: Jan 10, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_CURRENTS_AND_HEATING_H_
#define INCLUDE_CURRENTS_AND_HEATING_H_

/* ---------------------------------------------------------------------
 *
 * A basic current and heat equation solver with protrusion-like geometry
 *
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.templates.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>
#include <deal.II/lac/sparse_direct.h>	// UMFpack

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <set>


#include "laplace.h"

namespace Emitter {

using namespace dealii;

class CurrentsAndHeating {
public:
	CurrentsAndHeating();

	void run();

private:
	void make_grid();
	void setup_system();
	void assemble_system();
	void assemble_system_newton();
	void solve();
	void output_results(const unsigned int iteration) const;

	double em_current(const Point<2> &p, double t);

	Laplace laplace_problem;

	Triangulation<2> triangulation;
	FESystem<2> fe;
	DoFHandler<2> dof_handler;

	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double> system_rhs;

	Vector<double> solution; // u_{k+1}
	Vector<double> previous_solution; // u_k

	enum BoundaryId {
		not_specified = 0, // all boundaries which will have natural BCs
		copper_boundary = 1,
		bulk_bottom = 2
	};

	static constexpr double radius = 1.0;
	static constexpr double height = 4.0;
};

}



#endif /* INCLUDE_CURRENTS_AND_HEATING_H_ */
