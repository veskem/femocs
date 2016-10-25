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
CurrentsAndHeating<dim>::CurrentsAndHeating(PhysicalQuantities pq_, Laplace<dim>* laplace_,
		   double tip_temp_prediction_, double tip_pot_prediction_) :
		tip_temp_prediction(tip_temp_prediction_),
		tip_pot_prediction(tip_pot_prediction_),
		fe (FE_Q<dim>(currents_degree), 1, 	// Finite element type (1) = linear, etc and number of components
			FE_Q<dim>(heating_degree), 1),	// (we have 2 variables: potential and T with 1 component each)
		dof_handler(triangulation),
		pq(pq_),
		laplace(laplace_) {}

template <int dim>
Triangulation<dim>* CurrentsAndHeating<dim>::getp_triangulation() {
	return &triangulation;
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

	newton_update.reinit(dof_handler.n_dofs());
	present_solution.reinit(dof_handler.n_dofs());
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
		for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; f++) {
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
		for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cop_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
				Point<dim> cop_face_center = cop_cell->face(f)->center();
				// Loop over vacuum side and find corresponding (cell, face) pair
				for (unsigned int i=0; i < vacuum_interface_indexes.size(); i++) {
					if (cop_face_center.distance(vacuum_interface_centers[i]) < eps) {
						std::pair< std::pair<unsigned, unsigned>,
								   std::pair<unsigned, unsigned> > pair;
						pair.first = std::pair<unsigned, unsigned>(cop_cell->index(), f);
						pair.second = std::pair<unsigned, unsigned>(vacuum_interface_indexes[i], vacuum_interface_face[i]);
						interface_map.insert(pair);
					}
				}
			}
		}
	}
	// ---------------------------------------------------------------------------------------------

}


template <int dim>
class InitialValues : public Function<dim> {
	double tip_pot, tip_temp, amb_temp, zmax, zmin;
public:
	InitialValues (double tip_pot, double tip_temp, double amb_temp, double zmax, double zmin)
		: Function<dim>(2), tip_pot(tip_pot), tip_temp(tip_temp), amb_temp(amb_temp), zmax(zmax), zmin(zmin) {}
	double value (const Point<dim> &p,
				  const unsigned int component = 0) const {
		double z = 0.0;
		if (dim == 2) z = p(1);
		else z = p(2);

		double x = (z-zmin)/(zmax-zmin);

		if (component == 0) return tip_pot*x*x*x*x;
		else return amb_temp+(tip_temp-amb_temp)*x*x*x*x;
	}
};

template <int dim>
void CurrentsAndHeating<dim>::set_initial_condition() {

	typename Triangulation<dim>::active_face_iterator face;
	double zmax = -1e16, zmin = 1e16;

	// Loop through the faces and find maximum and minimum values for coordinates
	for (face = triangulation.begin_face(); face != triangulation.end_face(); ++face) {
		if (face->at_boundary()) {
			double z = face->center()[dim-1];
			if (z > zmax) zmax = z;
			if (z < zmin) zmin = z;
		}
	}
	// Set initial values such that the dirichlet BCs are satisfied
	// and the tip has the  predicted potential and temperature values.
	VectorTools::interpolate(dof_handler,
			InitialValues<dim>(tip_pot_prediction, tip_temp_prediction, ambient_temperature, zmax, zmin),
			present_solution);
}


// Assembles the linear system for one Newton iteration
template <int dim>
void CurrentsAndHeating<dim>::assemble_system_newton(const bool first_iteration) {

	TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);
	timer.enter_section("Pre-assembly");

	QGauss<dim> quadrature_formula(std::max(currents_degree, heating_degree)+1);
	QGauss<dim-1> face_quadrature_formula(std::max(std::max(currents_degree,
			heating_degree), laplace->shape_degree)+1);

	FEValues<dim> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values);

	FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
			update_values | update_gradients | update_normal_vectors
					| update_quadrature_points | update_JxW_values);

	// For evaluating E field at the copper surface
	FEFaceValues<dim> vacuum_fe_face_values(laplace->fe, face_quadrature_formula,
				update_gradients | update_quadrature_points | update_JxW_values);


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
	std::vector<Tensor<1,dim>> prev_sol_face_temperature_gradients (n_face_q_points);

	// Shape function values and gradients (arrays for every cell DOF)
	std::vector<double>				potential_phi (dofs_per_cell);
	std::vector< Tensor<1, dim> >	potential_phi_grad (dofs_per_cell);
	std::vector<double> 			temperature_phi (dofs_per_cell);
	std::vector< Tensor<1, dim> >	temperature_phi_grad (dofs_per_cell);

	// Electric field values from laplace solver
	std::vector< Tensor<1, dim> > electric_field_values (n_face_q_points);
	// ----------------------------------------------------------------------------------------------

    const FEValuesExtractors::Scalar potential (0);
    const FEValuesExtractors::Scalar temperature (1);

    timer.exit_section();
    timer.enter_section("Loop header");

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {

		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		fe_values[potential].get_function_gradients(present_solution, prev_sol_potential_gradients);
		fe_values[temperature].get_function_values(present_solution, prev_sol_temperature_values);
		fe_values[temperature].get_function_gradients(present_solution, prev_sol_temperature_gradients);

		timer.exit_section();
		timer.enter_section("Matrix assembly 1");

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

			timer.exit_section();
			timer.enter_section("Matrix assembly 2");

			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j) {
					cell_matrix(i, j)
						+= ( - (potential_phi_grad[i] * sigma * potential_phi_grad[j])
							 - (potential_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
							 + (temperature_phi[i] * 2 * sigma * prev_pot_grad * potential_phi_grad[j])
							 + (temperature_phi[i] * dsigma * prev_pot_grad * prev_pot_grad * temperature_phi[j])
							 - (temperature_phi_grad[i] * dkappa * prev_temp_grad * temperature_phi[j])
							 - (temperature_phi_grad[i] * kappa * temperature_phi_grad[j])
						   ) * fe_values.JxW(q);
					cell_matrix(i, j)
						+= ( - (potential_phi_grad[i] * sigma * potential_phi_grad[j])
							 - (potential_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
							 + (temperature_phi[i] * 2 * sigma * prev_pot_grad * potential_phi_grad[j])
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
			timer.exit_section();
			timer.enter_section("Matrix assembly 1");
		}
		// ---------------------------------------------------------------------------------------------
		timer.exit_section();
		timer.enter_section("Rhs assembly");
		// ---------------------------------------------------------------------------------------------
		// Local right-hand side assembly
		// ---------------------------------------------------------------------------------------------
		// integration over the boundary (cell faces)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary()) {
				fe_face_values.reinit(cell, f);

				fe_face_values[potential].get_function_gradients(present_solution, prev_sol_face_potential_gradients);
				fe_face_values[temperature].get_function_values(present_solution, prev_sol_face_temperature_values);
				fe_face_values[temperature].get_function_gradients(present_solution, prev_sol_face_temperature_gradients);

				if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
					// ---------------------------------------------------------------------------------------------
					// Vacuum side stuff
					// find the corresponding vacuum side face to the copper side face
					std::pair<unsigned, unsigned> vac_cell_info =
							interface_map[std::pair<unsigned, unsigned>(cell->index(), f)];
					// Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
					typename DoFHandler<dim>::active_cell_iterator vac_cell(&(laplace->triangulation),
							0, vac_cell_info.first, &(laplace->dof_handler));

					vacuum_fe_face_values.reinit(vac_cell, vac_cell_info.second);
					vacuum_fe_face_values.get_function_gradients(laplace->solution, electric_field_values);
					// ---------------------------------------------------------------------------------------------

					// loop through the quadrature points
					for (unsigned int q = 0; q < n_face_q_points; ++q) {

						double prev_temp = prev_sol_face_temperature_values[q];
						const Tensor<1,dim> prev_pot_grad = prev_sol_face_potential_gradients[q];
						const Tensor<1,dim> prev_temp_grad = prev_sol_face_temperature_gradients[q];

						const Tensor<1,dim> normal_vector = fe_face_values.normal_vector(q);

						//double sigma = pq.sigma(prev_temp);
						//double kappa = pq.kappa(prev_temp);
						double dsigma = pq.dsigma(prev_temp);
						double dkappa = pq.dkappa(prev_temp);
						double e_field = electric_field_values[q].norm();
						double emission_current = pq.emission_current(e_field, prev_temp);
						// Nottingham heat flux in
						// (eV*A/nm^2) -> (eV*n*q_e/(s*nm^2)) -> (J*n/(s*nm^2)) -> (W/nm^2)
						double nottingham_flux = pq.nottingham_de(e_field, prev_temp)*emission_current;

						for (unsigned int k=0; k<dofs_per_cell; ++k) {
							potential_phi[k] = fe_face_values[potential].value(k, q);
							temperature_phi[k] = fe_face_values[temperature].value(k, q);
						}
						for (unsigned int i = 0; i < dofs_per_cell; ++i) {
							cell_rhs(i) += 	(- (potential_phi[i] * emission_current)
											 - (temperature_phi[i] * nottingham_flux)
											)* fe_face_values.JxW(q);

							/*
							std::cout << (potential_phi[i] * sigma * normal_vector * prev_pot_grad) << " "
									  << (potential_phi[i] * emission_current) << std::endl;
							*/
							for (unsigned int j = 0; j < dofs_per_cell; ++j) {
								cell_matrix(i, j) += ( (potential_phi[i] * normal_vector * dsigma
														 	  * prev_pot_grad * temperature_phi[j])
													  + (temperature_phi[i] * normal_vector * dkappa
															  * prev_temp_grad * temperature_phi[j])
													 ) * fe_face_values.JxW(q);
							}

						}
					}
				}
			}
		}
		// ---------------------------------------------------------------------------------------------

		timer.exit_section();
		timer.enter_section("Loop footer");

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
		timer.exit_section();
		timer.enter_section("Loop header");
	}

	timer.exit_section();
	timer.enter_section("Post assembly");

	// Setting Dirichlet boundary values //

	// 0 potential at the bulk bottom boundary
	std::map<types::global_dof_index, double> current_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom, ZeroFunction<dim>(2),
			current_dirichlet, fe.component_mask(potential));

	MatrixTools::apply_boundary_values(current_dirichlet, system_matrix,
			newton_update, system_rhs);


	// 300K at bulk bottom if initial step, 0 otherwise
	std::map<types::global_dof_index, double> temperature_dirichlet;

	if (first_iteration) {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
				ConstantFunction<dim>(ambient_temperature, 2), temperature_dirichlet, fe.component_mask(temperature));
	} else {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
				ZeroFunction<dim>(2), temperature_dirichlet, fe.component_mask(temperature));
	}

	MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix,
			newton_update, system_rhs);

	timer.exit_section();
}


template <int dim>
void CurrentsAndHeating<dim>::solve() {
	// CG doesn't work as the matrix is not symmetric

	// GMRES
	/*
	SolverControl solver_control(400000, 1e-9);
	SolverGMRES<> solver_gmres(solver_control, SolverGMRES<>::AdditionalData(50));
	solver_gmres.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
	std::cout << "   " << solver_control.last_step() << " GMRES iterations needed to obtain convergence." << std::endl;
	*/


	// UMFPACK solver
	deallog << "Solving linear system with UMFPACK... " << std::endl;
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (newton_update, system_rhs);
}

template <int dim>
void CurrentsAndHeating<dim>::run() {

	std::cout << "/---------------------------------------------------------------/" << std::endl
		      << "CurrentsAndHeating run():" << std::endl;

	double temperature_tolerance = 1.0;

	Timer timer;

	setup_system();
	setup_mapping();

	// Sets the initial state based on a prediction.
	// Dirichlet BCs need to hold for this state
	// and 0 dirichlet BC should be applied for all Newton iterations
	set_initial_condition();

	std::cout << "    Setup and IC: " << timer.wall_time() << " s" << std::endl;

	// Newton iterations
	for (unsigned int iteration=0; iteration<10; ++iteration) {
		std::cout << "/--------------------------------/" << std::endl;
		std::cout << "Newton iteration " << iteration << std::endl;

		timer.restart();
		// reset the state of the linear system
		system_matrix.reinit(sparsity_pattern);
		system_rhs.reinit(dof_handler.n_dofs());
		std::cout << "    Reset state: " << timer.wall_time() << " s" << std::endl; timer.restart();

		timer.restart();
		//assemble_system_newton(iteration == 0);
		// No need to set dirichlet BC-s because they are set in set_initial_condition
		assemble_system_newton(false);

		std::cout << "    Assembly: " << timer.wall_time() << " s" << std::endl; timer.restart();

		solve();
		present_solution.add(1.0, newton_update); // alpha = 1.0

		std::cout << "    Solver: " << timer.wall_time() << " s" << std::endl; timer.restart();

		output_results(iteration);
		std::cout << "    output_results: " << timer.wall_time() << " s" << std::endl; timer.restart();

		std::cout << "    ||u_k-u_{k-1}||_L2 = " << newton_update.l2_norm() << std::endl;
		std::cout << "    ||u_k-u_{k-1}||_Linf = " << newton_update.linfty_norm() << std::endl;
		std::cout << "    Residual = " << system_rhs.l2_norm() << std::endl;

		if (newton_update.linfty_norm() < temperature_tolerance) {
			std::cout << "    Maximum temperature change less than tolerance: converged!" << std::endl;
			break;
		}
	}
	std::cout << "/---------------------------------------------------------------/" << std::endl;
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
				double e_field = -duh[i][0][d]; // gradient of the 0-th vector (i.e. potential)
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
	data_out.add_data_vector(present_solution, solution_names);
	data_out.add_data_vector(present_solution, current_post_processor);
	data_out.add_data_vector(present_solution, sigma_post_processor);
	data_out.build_patches();

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

template class CurrentsAndHeating<2> ;
template class CurrentsAndHeating<3> ;

} // namespace currents_heating

