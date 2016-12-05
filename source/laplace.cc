


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

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
#include <deal.II/lac/precondition.h>

#include "laplace.h"

namespace fch {
using namespace dealii;


// ----------------------------------------------------------------------------------------
// Class for outputting the resulting field distribution (calculated from potential distr.)

template <int dim>
class LaplacePostProcessor : public DataPostprocessorVector<dim> {
public:
	LaplacePostProcessor() : DataPostprocessorVector<dim>("Field", update_values | update_gradients) {}

	void
	compute_derived_quantities_scalar (	const std::vector<double>            	&/*uh*/,
										const std::vector<Tensor<1,dim> >  		&duh,
										const std::vector<Tensor<2,dim> >  		&/*dduh*/,
										const std::vector<Point<dim> >        	&/*normals*/,
										const std::vector<Point<dim> >          &/*evaluation_points*/,
										std::vector<Vector<double> >          	&computed_quantities) const {
		for (unsigned int i=0; i<computed_quantities.size(); i++) {
		    for (unsigned int d=0; d<dim; ++d)
		    	computed_quantities[i](d) = duh[i][d];

		}
	}
};
// ----------------------------------------------------------------------------


template<int dim>
Laplace<dim>::Laplace() :
		applied_efield(applied_efield_default), fe(shape_degree), dof_handler(triangulation) {
}

template<int dim>
Triangulation<dim>* Laplace<dim>::get_triangulation() {
	return &triangulation;
}

template <int dim>
DoFHandler<dim>* Laplace<dim>::get_dof_handler() {
	return &dof_handler;
}


template<int dim>
void Laplace<dim>::set_applied_efield(const double applied_field_) {
	applied_efield = applied_field_;
}


template<int dim>
void Laplace<dim>::import_mesh_from_file(const std::string file_name) {
	MeshPreparer<dim> mesh_preparer;

	mesh_preparer.import_mesh_from_file(&triangulation, file_name);
	mesh_preparer.mark_vacuum_boundary(&triangulation);
}

template<int dim>
bool Laplace<dim>::import_mesh_directly(std::vector<Point<dim> > vertices, std::vector<CellData<dim> > cells) {

	try {
		SubCellData subcelldata;
		// Do some clean-up on vertices...
		GridTools::delete_unused_vertices(vertices, cells, subcelldata);
		// ... and on cells
		GridReordering<dim, dim>::invert_all_cells_of_negative_grid(vertices, cells);

		triangulation.create_triangulation_compatibility(vertices, cells, SubCellData());
	} catch (exception &exc) {
		return false;
	}

	MeshPreparer<dim> mesh_preparer;
	mesh_preparer.mark_vacuum_boundary(&triangulation);

	return true;
}

template<int dim>
void Laplace<dim>::output_mesh(const std::string file_name) {
	MeshPreparer<dim> mesh_preparer;
	mesh_preparer.output_mesh(&triangulation, file_name);
}

template<int dim>
double Laplace<dim>::probe_efield(const Point<dim> &p) const {
	return VectorTools::point_gradient (dof_handler, solution, p).norm();
}

template<int dim>
std::vector<double> Laplace<dim>::get_potential(const std::vector<int> &cell_indexes,
												const std::vector<int> &vert_indexes) {

	// Initialise potentials with a value that is immediately visible if it's not changed to proper one
	std::vector<double> potentials(cell_indexes.size(), 1e15);

	for (unsigned i = 0; i < cell_indexes.size(); i++) {
		// Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
		// NB: only works without refinement !!!
		typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation,
				0, cell_indexes[i], &dof_handler);

		double potential = solution[dof_cell->vertex_dof_index(vert_indexes[i], 0)];
		potentials[i] = potential;
	}
    return potentials;
}

template<int dim>
std::vector<Tensor<1, dim> > Laplace<dim>::get_efield(const std::vector<int> &cell_indexes,
													 const std::vector<int> &vert_indexes) {

	QGauss<dim> quadrature_formula(quadrature_degree);
	FEValues<dim> fe_values(fe, quadrature_formula, update_gradients);

	std::vector< Tensor<1, dim> > solution_gradients (quadrature_formula.size());

	std::vector<Tensor<1, dim> > fields(cell_indexes.size());

	for (unsigned i = 0; i < cell_indexes.size(); i++) {
		// Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
		// NB: only works without refinement !!!
		typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation,
				0, cell_indexes[i], &dof_handler);

		fe_values.reinit(dof_cell);
		fe_values.get_function_gradients(solution, solution_gradients);
		Tensor<1, dim> field = -1.0 * solution_gradients.at(vert_indexes[i]);

		fields[i] = field;
	}
    return fields;
}

template<int dim>
void Laplace<dim>::setup_system() {
	dof_handler.distribute_dofs(fe);

	//std::cout << "    Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);

	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
void Laplace<dim>::assemble_system() {
	QGauss<dim> quadrature_formula(quadrature_degree);
	QGauss<dim-1> face_quadrature_formula(quadrature_degree);

	FEValues<dim> fe_values(fe, quadrature_formula,
			update_gradients | update_quadrature_points | update_JxW_values);

	FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
				update_values | update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

	for (; cell != endc; ++cell) {
		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q = 0; q < n_q_points; ++q) {
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
					cell_matrix(i, j) += (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
							* fe_values.JxW(q));

				//cell_rhs(i) += (fe_values.shape_value(i, q_index) * right_hand_side * fe_values.JxW(q_index));
			}
		}

		// Neumann boundary condition at faces with BoundaryId::vacuum_top
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == BoundaryId::vacuum_top) {
				fe_face_values.reinit(cell, f);

				for (unsigned int q = 0; q < n_face_q_points; ++q) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						cell_rhs(i) += (fe_face_values.shape_value(i, q)
								* applied_efield * fe_face_values.JxW(q));
					}
				}
			}
		}

		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_surface,
			ZeroFunction<dim>(), boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
}

template<int dim>
void Laplace<dim>::solve() {
	SolverControl solver_control(2000, 1e-9);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

	//std::cout << "   " << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;
}

template<int dim>
void Laplace<dim>::output_results(const std::string filename) const {
	LaplacePostProcessor<dim> field_calculator; // needs to be before data_out
	DataOut<dim> data_out;

	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "potential_v");
	data_out.add_data_vector(solution, field_calculator);

	data_out.build_patches();

	try {
		std::ofstream output(filename);
		data_out.write_vtk(output);
	} catch (...) {
		std::cerr << "WARNING: Couldn't open " + filename << ". ";
		std::cerr << "Output is not saved." << std::endl;
	}
}

template<int dim>
void Laplace<dim>::run() {
	Timer timer;
	std::cout << "/---------------------------------------------------------------/" << std::endl;
	std::cout << "Laplace solver: " << std::endl;
	setup_system();
	std::cout << "    setup_system(): " << timer.wall_time() << " s" << std::endl; timer.restart();
	assemble_system();
	std::cout << "    assemble_system(): " << timer.wall_time() << " s" << std::endl; timer.restart();
	solve();
	std::cout << "    solve(): " << timer.wall_time() << " s" << std::endl; timer.restart();
	output_results();
	std::cout << "    output_results(): " << timer.wall_time() << " s" << std::endl; timer.restart();
	std::cout << "/---------------------------------------------------------------/" << std::endl;

}

template class Laplace<2> ;
template class Laplace<3> ;

} // namespace fch
