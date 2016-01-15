/*
 * laplace.h
 *
 *  Created on: Jan 09, 2016
 *      Author: kristjan
 */

#ifndef LAPLACE_H_
#define LAPLACE_H_

/* ---------------------------------------------------------------------
 *
 * A basic laplace solver with protrusion-like geometry
 *
 */

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.templates.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <set>

namespace Emitter {
using namespace dealii;

class Laplace {
public:
	Laplace();

	void run();

	double probe_field(const Point<2> &p) const;

private:
	void make_grid_old();
	void make_grid();
	void setup_system();
	void assemble_system();
	void solve();
	void output_results() const;

	Triangulation<2> triangulation;
	FE_Q<2> fe;
	DoFHandler<2> dof_handler;

	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution;
	Vector<double> system_rhs;

	enum BoundaryId {
		not_specified = 0, // all boundaries which will have natural BCs
		copper_boundary = 1,
		vacuum_top = 2
	};
};

}

#endif /* LAPLACE_H_ */
