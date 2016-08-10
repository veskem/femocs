/*
 * currents_and_heating.cc
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#include "currents_and_heating.h"


template <int dim>
CurrentsAndHeating<dim>::CurrentsAndHeating(PhysicalQuantities pq_, Laplace<dim>* laplace_) :
		fe (FE_Q<dim>(1), 1, 	// Finite element type (1) = linear, etc and number of components
			FE_Q<dim>(1), 1),	// (we have 2 variables: potential and T with 1 component each)
		dof_handler(triangulation),
		pq(pq_),
		laplace(laplace_) {
}

template <int dim>
Triangulation<dim>* CurrentsAndHeating<dim>::getp_triangulation() {
	return &triangulation;
}

template <int dim>
double CurrentsAndHeating<dim>::emission_at_point(const Point<dim> &/*p*/, double temperature) {
	//double e_field = laplace->probe_field(p);
	double e_field = 10.0;
	return pq.emission_current(e_field, temperature);
}


template <int dim>
void CurrentsAndHeating<dim>::setup_system() {
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


// Assembles the linear system for one Newton iteration
template <int dim>
void CurrentsAndHeating<dim>::assemble_system_newton(const bool first_iteration) {

	Timer timer;
	timer.start();

	QGauss<dim> quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);

	FEValues<dim> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values);

	FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
			update_values | update_gradients | update_normal_vectors
					| update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	// The previous solution values in the cell quadrature points
	std::vector<Tensor<1,dim>> prev_sol_potential_gradients (n_q_points);
	std::vector<double> prev_sol_temperature_values (n_q_points);
	std::vector<Tensor<1,dim>> prev_sol_temperature_gradients (n_q_points);

	// The previous solution values in the face quadrature points
	std::vector<Tensor<1,dim>> prev_sol_face_potential_gradients (n_face_q_points);
	std::vector<double> prev_sol_face_temperature_values (n_face_q_points);

    const FEValuesExtractors::Scalar potential (0);
    const FEValuesExtractors::Scalar temperature (1);

    timer.stop();
    std::cout << "Initialization: " << timer() << std::endl;

    double total_matrix_time = 0.0;
    double total_faces_time = 0.0;
    double total_add_global_time = 0.0;

    int field_probe_count = 0;

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {

		timer.restart();
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		fe_values[potential].get_function_gradients(previous_solution, prev_sol_potential_gradients);
		fe_values[temperature].get_function_values(previous_solution, prev_sol_temperature_values);
		fe_values[temperature].get_function_gradients(previous_solution, prev_sol_temperature_gradients);

		for (unsigned int q = 0; q < n_q_points; ++q) {
			double prev_temp = prev_sol_temperature_values[q];
			const Tensor<1,dim> prev_pot_grad = prev_sol_potential_gradients[q];
			const Tensor<1,dim> prev_temp_grad = prev_sol_temperature_gradients[q];

			double sigma = pq.sigma(prev_temp);
			double dsigma = pq.dsigma(prev_temp);
			double kappa = pq.kappa(prev_temp);
			double dkappa = pq.dkappa(prev_temp);

			for (unsigned int i = 0; i < dofs_per_cell; ++i) {

				// these should be evaluated outside the loops (into an array)
				const Tensor<1,dim> grad_phi_i_p	= fe_values[potential].gradient(i, q);
				const double phi_i_t 				= fe_values[temperature].value(i, q);
				const Tensor<1,dim> grad_phi_i_t	= fe_values[temperature].gradient(i, q);

				for (unsigned int j = 0; j < dofs_per_cell; ++j) {

					const Tensor<1,dim> grad_phi_j_p	= fe_values[potential].gradient(j, q);
					const double phi_j_t 				= fe_values[temperature].value(j, q);
					const Tensor<1,dim> grad_phi_j_t	= fe_values[temperature].gradient(j, q);

					cell_matrix(i, j) += (	- (grad_phi_i_p * sigma * grad_phi_j_p)
											- (grad_phi_i_p * dsigma * prev_pot_grad * phi_j_t)
											+ (phi_i_t * sigma * 2 * prev_pot_grad * grad_phi_j_p)
											+ (phi_i_t * dsigma * prev_pot_grad * prev_pot_grad * phi_j_t)
											- (grad_phi_i_t * dkappa * prev_temp_grad * phi_j_t)
											- (grad_phi_i_t * kappa * grad_phi_j_t)
										 ) * fe_values.JxW(q);
				}
				cell_rhs(i) += (	grad_phi_i_p * sigma * prev_pot_grad
									- phi_i_t * sigma * prev_pot_grad * prev_pot_grad
									+ grad_phi_i_t * kappa * prev_temp_grad
								) * fe_values.JxW(q);
			}
		}


		timer.stop();
		total_matrix_time += timer();
		timer.restart();

		// integration over the boundary (cell faces)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary()) {
				fe_face_values.reinit(cell, f);

				fe_face_values[potential].get_function_gradients(previous_solution, prev_sol_face_potential_gradients);
				fe_face_values[temperature].get_function_values(previous_solution, prev_sol_face_temperature_values);

				if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
					// loop thought the quadrature points
					for (unsigned int q = 0; q < n_face_q_points; ++q) {

						double prev_temp = prev_sol_face_temperature_values[q];
						const Tensor<1,dim> prev_pot_grad = prev_sol_face_potential_gradients[q];

						double dsigma = pq.dsigma(prev_temp);
						field_probe_count++;
						double emission_current = emission_at_point(fe_face_values.quadrature_point(q), prev_temp);

						for (unsigned int i = 0; i < dofs_per_cell; ++i) {
							const double phi_i_p = fe_face_values[potential].value(i, q);
							// neumann BC for the current, one from LHS and one from RHS
							cell_rhs(i) += 	(- (phi_i_p * emission_current)
											)* fe_face_values.JxW(q);

							for (unsigned int j = 0; j < dofs_per_cell; ++j) {
								const double phi_j_t = fe_face_values[temperature].value(j, q);

								cell_matrix(i, j) += (	phi_i_p * fe_face_values.normal_vector(q)
														* prev_pot_grad * dsigma * phi_j_t
													 ) * fe_face_values.JxW(q);
							}

						}
					}
				}
			}

		timer.stop();
		total_faces_time += timer();
		timer.restart();

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}

		timer.stop();
		total_add_global_time += timer();
	}
	timer.restart();
	// Setting Dirichlet boundary values //

	// 0 potential at the bulk bottom boundary
	std::map<types::global_dof_index, double> current_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom, ZeroFunction<dim>(2),
			current_dirichlet, fe.component_mask(potential));

	MatrixTools::apply_boundary_values(current_dirichlet, system_matrix,
			solution, system_rhs);


	// 300K at bulk bottom if initial step, 0 otherwise
	std::map<types::global_dof_index, double> temperature_dirichlet;

	if (first_iteration) {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
				ConstantFunction<dim>(300.0, 2), temperature_dirichlet, fe.component_mask(temperature));
	} else {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
				ZeroFunction<dim>(2), temperature_dirichlet, fe.component_mask(temperature));
	}

	MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix,
			solution, system_rhs);

	timer.stop();

	std::cout << "Total matrix: " << total_matrix_time << "s" << std::endl;
	std::cout << "Total faces: " << total_faces_time << "s" << std::endl;
	std::cout << "Total global addition: " << total_add_global_time << "s" << std::endl;
	std::cout << "End stuff: " << timer() << "s" << std::endl;
	std::cout << "Probe count: " << field_probe_count << std::endl;
}


template <int dim>
void CurrentsAndHeating<dim>::solve() {
	// CG doesn't work as the matrix is not symmetric

	// UMFPACK solver
	deallog << "Solving linear system with UMFPACK... ";
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (solution, system_rhs);
}

template <int dim>
void CurrentsAndHeating<dim>::run() {

	setup_system();

	// Newton iterations
	for (unsigned int iteration=0; iteration<5; ++iteration) {
		std::cout << "    Newton iteration " << iteration << std::endl;

		// reset the state of the linear system
		system_matrix.reinit(sparsity_pattern);
		system_rhs.reinit(dof_handler.n_dofs());

		Timer timer;
		timer.start ();
		assemble_system_newton(iteration == 0);
		timer.stop ();

		deallog << "system assembled in "
				<< timer ()
				<< "s"
				<< std::endl;

		timer.restart();
		solve();
		// u_{k+1} = \alpha * \delta
		if (iteration!=0)
			solution *= 1.0; // alpha
		solution.add(1.0, previous_solution);
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


// ----------------------------------------------------------------------------------------
// Class for outputting the current density distribution (calculated from potential distr.)
template <int dim>
class CurrentPostProcessor : public DataPostprocessorVector<dim> {
	PhysicalQuantities pq;
public:
	CurrentPostProcessor(PhysicalQuantities pq_)
		: DataPostprocessorVector<dim>("current_density", update_values | update_gradients),
		  pq(pq_){}

	void
	compute_derived_quantities_vector (	const std::vector< Vector< double > > &  					uh,
											const std::vector< std::vector< Tensor< 1, dim > > > &  duh,
											const std::vector< std::vector< Tensor< 2, dim > > > &  /*dduh*/,
											const std::vector< Point< dim > > &  					/*normals*/,
											const std::vector< Point< dim > > &  					/*evaluation_points*/,
											std::vector< Vector< double > > &  						computed_quantities
										  ) const {
		for (unsigned int i=0; i<computed_quantities.size(); i++) {
			double t = uh[i][1]; // temperature
			double sigma = pq.sigma(t);
			for (unsigned int d=0; d<dim; ++d) {
				double e_field = duh[i][0][d]; // gradient of the 0-st vector (i.e. potential)
				computed_quantities[i](d) = sigma*e_field;
			}
		}
	}
};
// Class for outputting the electrical conductivity distribution
template <int dim>
class SigmaPostProcessor : public DataPostprocessorScalar<dim> {
	PhysicalQuantities pq;
public:
	SigmaPostProcessor(PhysicalQuantities pq_)
		: DataPostprocessorScalar<dim>("sigma", update_values),
		  pq(pq_){}

	void
	compute_derived_quantities_vector (	const std::vector< Vector< double > > &  					uh,
											const std::vector< std::vector< Tensor< 1, dim > > > &  /*duh*/,
											const std::vector< std::vector< Tensor< 2, dim > > > &  /*dduh*/,
											const std::vector< Point< dim > > &  					/*normals*/,
											const std::vector< Point< dim > > &  					/*evaluation_points*/,
											std::vector< Vector< double > > &  						computed_quantities
										  ) const {
		for (unsigned int i=0; i<computed_quantities.size(); i++) {
			double t = uh[i][1]; // temperature
			computed_quantities[i](0) = pq.sigma(t);
		}
	}
};
// ----------------------------------------------------------------------------------------

template <int dim>
void CurrentsAndHeating<dim>::output_results(const unsigned int iteration) const {

	std::string filename = "solution-" + Utilities::int_to_string(iteration) + ".vtk";

	std::vector<std::string> solution_names {"potential", "temperature"};

	CurrentPostProcessor<dim> current_post_processor(pq); // needs to be before data_out
	SigmaPostProcessor<dim> sigma_post_processor(pq); // needs to be before data_out

	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names);
	data_out.add_data_vector(solution, current_post_processor);
	data_out.add_data_vector(solution, sigma_post_processor);
	data_out.build_patches();

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

template class CurrentsAndHeating<2> ;
template class CurrentsAndHeating<3> ;
