

#include "currents_and_heating.h"
#include "physics.h"
#include "utility.h"

namespace Emitter {

CurrentsAndHeating::CurrentsAndHeating() :
		fe (FE_Q<2>(1), 1, 	// Finite element type (Q = usual ?) (1)
			FE_Q<2>(1), 1),	// additionally, for how many variables
		dof_handler(triangulation) {
}

// ----------------------------------------------------------------------------------------
// Function classes
// ----------------------------------------------------------------------------------------

// current density calculator
class CurrentPostProcessor : public DataPostprocessorVector<2> {
public:
	CurrentPostProcessor();

	virtual void
	compute_derived_quantities_vector (	const std::vector< Vector< double > > &  	uh,
										const std::vector< std::vector< Tensor< 1, 2 > > > &  	duh,
										const std::vector< std::vector< Tensor< 2, 2 > > > &  	dduh,
										const std::vector< Point< 2 > > &  	normals,
										const std::vector< Point< 2 > > &  	evaluation_points,
										std::vector< Vector< double > > &  	computed_quantities
									  ) const;
};

CurrentPostProcessor::CurrentPostProcessor() : DataPostprocessorVector<2>("current_density", update_values | update_gradients) {}

void CurrentPostProcessor::
compute_derived_quantities_vector (	const std::vector< Vector< double > > &  	uh,
										const std::vector< std::vector< Tensor< 1, 2 > > > &  	duh,
										const std::vector< std::vector< Tensor< 2, 2 > > > &  	/*dduh*/,
										const std::vector< Point< 2 > > &  	/*normals*/,
										const std::vector< Point< 2 > > &  	/*evaluation_points*/,
										std::vector< Vector< double > > &  	computed_quantities
									  ) const {
	for (unsigned int i=0; i<computed_quantities.size(); i++) {
		double t = uh[i][1]; // temperature
		for (unsigned int d=0; d<2; ++d) {
			double e_field = duh[i][0][d]; // gradient of the 0-st vector (i.e. potential)
			computed_quantities[i](d) = el_conductivity(t)*e_field;
		}
	}
}


template <int dim>
class AmbientTemperature : public Function<dim> {
public:
	AmbientTemperature() : Function<dim>(2) {} // The number of components needs to be passed to Function, probably...
	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
};
template <int dim>
double AmbientTemperature<dim>::value (const Point<dim> &/*p */, const unsigned int /*component */) const {
	return 300;
}

double CurrentsAndHeating::em_current(const Point<2> &p, double temperature) {
	double e_field = laplace_problem.probe_field(p);
	return emission_current(e_field, temperature);
}

// ----------------------------------------------------------------------------------------


void CurrentsAndHeating::make_grid() {

	double r = radius;	//radius of the tip
	double h = height;	//height of the tip

	// the bulk dimensions should be given in terms of r
	// to make the mesh points line up
	unsigned int bulk_r = 9; // this should be odd
	unsigned int bulk_h = 6; // this should be even

	double eps = 1e-6;

	unsigned char tip_manifold = 10;

	// --------------------------------------
	// Create a custom half-circle
	// --------------------------------------

	const Point<2> vertices_1[] = {
			Point<2>(-r, 0),
			Point<2>(-r/3, 0),
			Point<2>( r/3, 0),
			Point<2>( r, 0),
			Point<2>(-r/3, r/3*std::sqrt(3)),
			Point<2>( r/3, r/3*std::sqrt(3)),
			Point<2>(-r/2, r/2*std::sqrt(3)),
			Point<2>( r/2, r/2*std::sqrt(3))
	};
	const unsigned int n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
	const std::vector<Point<2> > vertices(&vertices_1[0],&vertices_1[n_vertices]);
	const int cell_vertices[][GeometryInfo<2>::vertices_per_cell] = {
			{0,1,6,4},
			{1,2,4,5},
			{3,7,2,5},
			{6,4,7,5}
	};
	const unsigned int n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);
	std::vector<CellData<2> > cells(n_cells, CellData<2>());
	for (unsigned int i = 0; i < n_cells; ++i) {
		for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
			cells[i].vertices[j] = cell_vertices[i][j];
		cells[i].material_id = 0;
	}
	triangulation.create_triangulation(vertices, cells, SubCellData());

	// --------------------------------------
	// Create the shaft of the protrusion
	// --------------------------------------

	Triangulation<2> mesh_shaft;

	unsigned int shaft_rep = (unsigned int)(1.5*h/r+1);
	std::vector<unsigned int> rep {3, shaft_rep};

	GridGenerator::subdivided_hyper_rectangle(mesh_shaft, rep,
			Point<2>(-r, -h), Point<2>(r, 0));
	GridGenerator::merge_triangulations(triangulation, mesh_shaft, triangulation);

	// --------------------------------------
	// Create the bulk copper
	// --------------------------------------

	Triangulation<2> mesh_bulk1;

	std::vector<unsigned int> rep1 {3*bulk_r,3*bulk_h/2};
	GridGenerator::subdivided_hyper_rectangle(mesh_bulk1, rep1,
				Point<2>(-(bulk_r*r), -h-(bulk_h*r)), Point<2>((bulk_r*r), -h));
	GridGenerator::merge_triangulations(triangulation, mesh_bulk1, triangulation);

	// --------------------------------------
	// Define manifolds and boundaries
	// --------------------------------------

	Triangulation<2>::active_cell_iterator cell;
	for (cell = triangulation.begin_active(); cell != triangulation.end(); ++cell) {
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
			// The face is on the spherical tip
			if (cell->face(f)->at_boundary()) {
				Point<2> center = cell->face(f)->center();

				// boundaries on the tip
				if (center[1] + eps > 0) {
					cell->face(f)->set_all_manifold_ids(tip_manifold);
				}
				// copper surface
				if (center[1] + eps > -h) {
					cell->face(f)->set_all_boundary_ids(copper_boundary);
				}
				// Bulk bottom
				if (center[1] - eps < -h-(bulk_h*r)) {
					cell->face(f)->set_all_boundary_ids(bulk_bottom);
				}

			}
		}
	}

	static const SphericalManifold<2> boundary;
	triangulation.set_manifold(tip_manifold, boundary);
	triangulation.refine_global(3);

	output_mesh(triangulation, "mesh.eps");
	triangulation.set_manifold(tip_manifold); // reset the spherical manifold
}

void CurrentsAndHeating::setup_system() {
	dof_handler.distribute_dofs(fe);
	std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
			<< std::endl;

	DoFRenumbering::component_wise (dof_handler);

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);

	solution.reinit(dof_handler.n_dofs());
	previous_solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

void CurrentsAndHeating::assemble_system() {
	QGauss<2> quadrature_formula(2);
	QGauss<1> face_quadrature_formula(2);

	FEValues<2> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values);

	FEFaceValues<2> fe_face_values(fe, face_quadrature_formula,
			update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	std::vector<Tensor<1,2>> prev_sol_potential_gradients (n_q_points);
	std::vector<double> prev_sol_temperature_values (n_q_points);

    const FEValuesExtractors::Scalar potential (0);
    const FEValuesExtractors::Scalar temperature (1);


	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		//electric_conductivity.value_list(fe_values.get_quadrature_points(), el_conductivity_values);

		fe_values[potential].get_function_gradients(previous_solution, prev_sol_potential_gradients);
		fe_values[temperature].get_function_values(previous_solution, prev_sol_temperature_values);
		// Better to calculate the conductivities once here...

		for (unsigned int q = 0; q < n_q_points; ++q) {
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {

				const Tensor<1,2> grad_phi_i_p	= fe_values[potential].gradient(i, q);
				const double phi_i_t 		= fe_values[temperature].value(i, q);
				const Tensor<1,2> grad_phi_i_t	= fe_values[temperature].gradient(i, q);

				for (unsigned int j = 0; j < dofs_per_cell; ++j) {

					const Tensor<1,2> grad_phi_j_p	= fe_values[potential].gradient(j, q);
					const Tensor<1,2> grad_phi_j_t	= fe_values[temperature].gradient(j, q);

					cell_matrix(i, j) += (	- (grad_phi_i_p * el_conductivity(prev_sol_temperature_values[q]) * grad_phi_j_p)
											+ (phi_i_t * el_conductivity(prev_sol_temperature_values[q])
												* prev_sol_potential_gradients[q] * grad_phi_j_p)
											- (grad_phi_i_t * th_conductivity(prev_sol_temperature_values[q]) * grad_phi_j_t)
										 ) * fe_values.JxW(q);

				}

				// RHS is 0 in our case...
				cell_rhs(i) += ( 0 * fe_values.JxW(q));
			}
		}

		// Neumann boundary conditions...
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary()
					&& cell->face(f)->boundary_id() == copper_boundary) {
				fe_face_values.reinit(cell, f);
				for (unsigned int q_index = 0; q_index < n_face_q_points;
						++q_index) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						// Current density BC on copper boundary
						cell_rhs(i) +=  - fe_face_values[potential].value(i, q_index)
												* em_current(fe_face_values.quadrature_point(q_index), prev_sol_temperature_values[q_index])
												* fe_face_values.JxW(q_index);
					}
				}
			}

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	// Setting Dirichlet boundary values //

	// 0 potential at the bulk bottom boundary
	std::map<types::global_dof_index, double> current_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, bulk_bottom, ZeroFunction<2>(2),
			current_dirichlet, fe.component_mask(potential));

	MatrixTools::apply_boundary_values(current_dirichlet, system_matrix,
			solution, system_rhs);

	// 300K at bulk bottom
	std::map<types::global_dof_index, double> temperature_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, bulk_bottom, AmbientTemperature<2>(),
			temperature_dirichlet, fe.component_mask(temperature));

	MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix,
			solution, system_rhs);
}

void CurrentsAndHeating::assemble_system_newton() {
	QGauss<2> quadrature_formula(2);
	QGauss<1> face_quadrature_formula(2);

	FEValues<2> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values);

	FEFaceValues<2> fe_face_values(fe, face_quadrature_formula,
			update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	std::vector<Tensor<1,2>> prev_sol_potential_gradients (n_q_points);
	std::vector<double> prev_sol_temperature_values (n_q_points);

    const FEValuesExtractors::Scalar potential (0);
    const FEValuesExtractors::Scalar temperature (1);


	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		//electric_conductivity.value_list(fe_values.get_quadrature_points(), el_conductivity_values);

		fe_values[potential].get_function_gradients(previous_solution, prev_sol_potential_gradients);
		fe_values[temperature].get_function_values(previous_solution, prev_sol_temperature_values);
		// Better to calculate the conductivities once here...

		for (unsigned int q = 0; q < n_q_points; ++q) {
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {

				const Tensor<1,2> grad_phi_i_p	= fe_values[potential].gradient(i, q);
				const double phi_i_t 		= fe_values[temperature].value(i, q);
				const Tensor<1,2> grad_phi_i_t	= fe_values[temperature].gradient(i, q);

				for (unsigned int j = 0; j < dofs_per_cell; ++j) {

					const Tensor<1,2> grad_phi_j_p	= fe_values[potential].gradient(j, q);
					const Tensor<1,2> grad_phi_j_t	= fe_values[temperature].gradient(j, q);

					cell_matrix(i, j) += (	- (grad_phi_i_p * el_conductivity(prev_sol_temperature_values[q]) * grad_phi_j_p)
											+ (phi_i_t * el_conductivity(prev_sol_temperature_values[q])
												* prev_sol_potential_gradients[q] * grad_phi_j_p)
											- (grad_phi_i_t * th_conductivity(prev_sol_temperature_values[q]) * grad_phi_j_t)
										 ) * fe_values.JxW(q);

				}

				// RHS is 0 in our case...
				cell_rhs(i) += ( 0 * fe_values.JxW(q));
			}
		}

		// Neumann boundary conditions...
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary()
					&& cell->face(f)->boundary_id() == copper_boundary) {
				fe_face_values.reinit(cell, f);
				for (unsigned int q_index = 0; q_index < n_face_q_points;
						++q_index) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						// Current density BC on copper boundary
						cell_rhs(i) +=  - fe_face_values[potential].value(i, q_index)
												* em_current(fe_face_values.quadrature_point(q_index), prev_sol_temperature_values[q_index])
												* fe_face_values.JxW(q_index);
					}
				}
			}

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	// Setting Dirichlet boundary values //

	// 0 potential at the bulk bottom boundary
	std::map<types::global_dof_index, double> current_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, bulk_bottom, ZeroFunction<2>(2),
			current_dirichlet, fe.component_mask(potential));

	MatrixTools::apply_boundary_values(current_dirichlet, system_matrix,
			solution, system_rhs);

	// 300K at bulk bottom
	std::map<types::global_dof_index, double> temperature_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, bulk_bottom, AmbientTemperature<2>(),
			temperature_dirichlet, fe.component_mask(temperature));

	MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix,
			solution, system_rhs);

}

void CurrentsAndHeating::solve() {

	/*
	// CG doesn't seem to work...
	deallog << "Solving linear system with CG... ";
	SolverControl solver_control(2000, 1e-12);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
	*/


	// UMFPACK solver
	deallog << "Solving linear system with UMFPACK... ";
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (solution, system_rhs);


}

void CurrentsAndHeating::output_results(const unsigned int iteration) const {

	std::string filename = "solution-" + Utilities::int_to_string(iteration) + ".vtk";

	std::vector<std::string> solution_names {"potential", "temperature"};

	CurrentPostProcessor field_calculator; // needs to be before data_out

	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names);
	data_out.add_data_vector(solution, field_calculator);
	data_out.build_patches();

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

void CurrentsAndHeating::run() {


	// calculate the field
	laplace_problem.run();

	make_grid();
	setup_system();

	// 10 Picard iterations
	for (unsigned int iteration=0; iteration<10; ++iteration) {

		std::cout << "    Iteration " << iteration << std::endl;

		Timer timer;
		timer.start ();
		assemble_system();
		timer.stop ();

		deallog << "system assembled in "
				<< timer ()
				<< "s"
				<< std::endl;

		timer.restart();
		solve();
		timer.stop ();
		deallog << "solver finished in "
				<< timer ()
				<< "s"
				<< std::endl;

		output_results(iteration);

		Vector<double> difference(dof_handler.n_dofs());
		difference = solution;
		difference -= previous_solution;
		std::cout << "    ||u_k-u_{k-1}|| = " << difference.l2_norm() << std::endl;

		previous_solution = solution;
	}

}

} // end namespace

