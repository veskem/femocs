/*
 * currents_and_heating.cc
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#include "currents_and_heating.h"

namespace currents_heating {
using namespace dealii;
using namespace laplace;

template <int dim>
CurrentsAndHeating<dim>::CurrentsAndHeating(PhysicalQuantities pq_, Laplace<dim>* laplace_) :
		fe (FE_Q<dim>(currents_degree), 1, 	// Finite element type (1) = linear, etc and number of components
			FE_Q<dim>(heating_degree), 1),	// (we have 2 variables: potential and T with 1 component each)
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

template <int dim>
void CurrentsAndHeating<dim>::setup_mapping() {

	double eps = 1e-12;

	// ---------------------------------------------------------------------------------------------
	// Loop over vacuum interface cells

	std::vector<unsigned int> vacuum_interface_indexes;
	std::vector<unsigned int> vacuum_interface_face;
	std::vector< Point<dim> > vacuum_interface_centers;

	typename DoFHandler<dim>::active_cell_iterator
		vac_cell = laplace->dof_handler.begin_active(),
		vac_endc = laplace->dof_handler.end();
	for (; vac_cell != vac_endc; ++vac_cell) {
		for (char f=0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (vac_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
				vacuum_interface_indexes.push_back(vac_cell->index());
				vacuum_interface_face.push_back(f);
				vacuum_interface_centers.push_back(vac_cell->face(f)->center());
			}
		}
	}
	// ---------------------------------------------------------------------------------------------

	// ---------------------------------------------------------------------------------------------
	// Loop over copper interface cells

	typename DoFHandler<dim>::active_cell_iterator
		cop_cell = dof_handler.begin_active(),
		cop_endc = dof_handler.end();

	for (; cop_cell != cop_endc; ++cop_cell) {
		for (char f=0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cop_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
				Point<dim> cop_face_center = cop_cell->face(f)->center();
				// Loop over vacuum side and find corresponding (cell, face) pair
				for (unsigned int i=0; i < vacuum_interface_indexes.size(); i++) {
					if (cop_face_center.distance(vacuum_interface_centers[i]) < eps) {
						std::pair< std::pair<unsigned int, char>, std::pair<unsigned int, char> > pair;
						pair.first = std::pair<unsigned int, char>(cop_cell->index(), f);
						pair.second = std::pair<unsigned int, char>(vacuum_interface_indexes[i], vacuum_interface_face[i]);
						interface_map.insert(pair);
					}
				}
			}
		}
	}
	// ---------------------------------------------------------------------------------------------

}


// Assembles the linear system for one Newton iteration
template <int dim>
void CurrentsAndHeating<dim>::assemble_system_newton(const bool first_iteration) {

	QGauss<dim> quadrature_formula(std::max(currents_degree, heating_degree)+2);
	QGauss<dim-1> face_quadrature_formula(std::max(currents_degree, heating_degree)+2);

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

	// ---------------------------------------------------------------------------------------------
	// Declaring variables once at the start (for performance reasons)

	// The previous solution values in the cell quadrature points
	std::vector<Tensor<1,dim>> prev_sol_potential_gradients (n_q_points);
	std::vector<double> prev_sol_temperature_values (n_q_points);
	std::vector<Tensor<1,dim>> prev_sol_temperature_gradients (n_q_points);

	// The previous solution values in the face quadrature points
	std::vector<Tensor<1,dim>> prev_sol_face_potential_gradients (n_face_q_points);
	std::vector<double> prev_sol_face_temperature_values (n_face_q_points);

	// Shape function values and gradients (arrays for every cell DOF)
	std::vector<double>				potential_phi (dofs_per_cell);
	std::vector< Tensor<1, dim> >	potential_phi_grad (dofs_per_cell);
	std::vector<double> 			temperature_phi (dofs_per_cell);
	std::vector< Tensor<1, dim> >	temperature_phi_grad (dofs_per_cell);
	// ----------------------------------------------------------------------------------------------

    const FEValuesExtractors::Scalar potential (0);
    const FEValuesExtractors::Scalar temperature (1);

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {

		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		fe_values[potential].get_function_gradients(previous_solution, prev_sol_potential_gradients);
		fe_values[temperature].get_function_values(previous_solution, prev_sol_temperature_values);
		fe_values[temperature].get_function_gradients(previous_solution, prev_sol_temperature_gradients);

		// ---------------------------------------------------------------------------------------------
		// Local matrix assembly
		// ---------------------------------------------------------------------------------------------
		for (unsigned int q = 0; q < n_q_points; ++q) {
			double prev_temp = prev_sol_temperature_values[q];
			const Tensor<1,dim> prev_pot_grad = prev_sol_potential_gradients[q];
			const Tensor<1,dim> prev_temp_grad = prev_sol_temperature_gradients[q];

			double sigma = pq.sigma(prev_temp);
			double dsigma = pq.dsigma(prev_temp);
			double kappa = pq.kappa(prev_temp);
			double dkappa = pq.dkappa(prev_temp);

			for (unsigned int k=0; k<dofs_per_cell; ++k) {
				potential_phi_grad[k] = fe_values[potential].gradient (k, q);
				temperature_phi[k] = fe_values[temperature].value (k, q);
				temperature_phi_grad[k] = fe_values[temperature].gradient (k, q);
			}
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j) {
					cell_matrix(i, j)
						+= ( - (potential_phi_grad[i] * sigma * potential_phi_grad[j])
							 - (potential_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
							 + (temperature_phi[i] * sigma * 2 * prev_pot_grad * potential_phi_grad[j])
							 + (temperature_phi[i] * dsigma * prev_pot_grad * prev_pot_grad * temperature_phi[j])
							 - (temperature_phi_grad[i] * dkappa * prev_temp_grad * temperature_phi[j])
							 - (temperature_phi_grad[i] * kappa * temperature_phi_grad[j])
						   ) * fe_values.JxW(q);
				}
				cell_rhs(i)
					+= (   potential_phi_grad[i] * sigma * prev_pot_grad
						 - temperature_phi[i] * sigma * prev_pot_grad * prev_pot_grad
						 + temperature_phi_grad[i] * kappa * prev_temp_grad
					   ) * fe_values.JxW(q);
			}
		}
		// ---------------------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------------------
		// Local right-hand side assembly
		// ---------------------------------------------------------------------------------------------
		// integration over the boundary (cell faces)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary()) {
				fe_face_values.reinit(cell, f);

				fe_face_values[potential].get_function_gradients(previous_solution, prev_sol_face_potential_gradients);
				fe_face_values[temperature].get_function_values(previous_solution, prev_sol_face_temperature_values);

				if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
					// ---------------------------------------------------------------------------------------------
					// Vacuum side stuff
					interface_map[std::pair<cell->index(), f>];

					// ---------------------------------------------------------------------------------------------

					// loop thought the quadrature points
					for (unsigned int q = 0; q < n_face_q_points; ++q) {

						double prev_temp = prev_sol_face_temperature_values[q];
						const Tensor<1,dim> prev_pot_grad = prev_sol_face_potential_gradients[q];

						double dsigma = pq.dsigma(prev_temp);
						double emission_current = emission_at_point(fe_face_values.quadrature_point(q), prev_temp);

						for (unsigned int k=0; k<dofs_per_cell; ++k) {
							potential_phi[k] = fe_values[potential].value (k, q);
							temperature_phi[k] = fe_values[temperature].value (k, q);
						}
						for (unsigned int i = 0; i < dofs_per_cell; ++i) {
							cell_rhs(i) += 	(- (potential_phi[i] * emission_current)
											)* fe_face_values.JxW(q);

							for (unsigned int j = 0; j < dofs_per_cell; ++j) {
								cell_matrix(i, j) += (	potential_phi[i] * fe_face_values.normal_vector(q)
														* prev_pot_grad * dsigma * temperature_phi[j]
													 ) * fe_face_values.JxW(q);
							}

						}
					}
				}
			}
		}
		// ---------------------------------------------------------------------------------------------


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

	setup_mapping();
/*
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
*/
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
				double e_field = duh[i][0][d]; // gradient of the 0-th vector (i.e. potential)
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

} // namespace currents_heating

