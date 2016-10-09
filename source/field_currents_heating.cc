/*
 * field_currents_heating.cc
 *
 *  Created on: Aug 3, 2016
 *      Author: kristjan
 */


#include "field_currents_heating.h"

namespace field_currents_heating {
using namespace dealii;


template <int dim>
FieldCurrentsHeating<dim>::FieldCurrentsHeating(PhysicalQuantities pq_)
		:
		vacuum_fe(FE_Q<dim>(field_degree), 1,
				  FE_Nothing<dim>(), 1,
				  FE_Nothing<dim>(), 1),
		copper_fe(FE_Nothing<dim>(), 1,
				  FE_Q<dim>(currents_degree), 1,
				  FE_Q<dim>(heating_degree), 1),
		dof_handler(triangulation),
		pq(pq_) {
	fe_collection.push_back(vacuum_fe);
	fe_collection.push_back(copper_fe);
}

template <int dim>
Triangulation<dim>* FieldCurrentsHeating<dim>::getp_triangulation() {
	return &triangulation;
}

template <int dim>
bool FieldCurrentsHeating<dim>::
cell_is_in_vacuum_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell) {
  return (cell->material_id() == MaterialId::vacuum_domain);
}

template <int dim>
bool FieldCurrentsHeating<dim>::
cell_is_in_copper_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell) {
  return (cell->material_id() == MaterialId::copper_domain);
}

template <int dim>
void FieldCurrentsHeating<dim>::set_active_fe_indices () {
  for (typename hp::DoFHandler<dim>::active_cell_iterator
	   cell = dof_handler.begin_active();
	   cell != dof_handler.end(); ++cell)
	{
	  if (cell_is_in_vacuum_domain(cell))
		cell->set_active_fe_index (0); // corresponds to fe_collection[0] = vacuum_fe
	  else if (cell_is_in_copper_domain(cell))
		cell->set_active_fe_index (1); // corresponds to fe_collection[1] = copper_fe
	  else
		Assert (false, ExcNotImplemented());
	}
}

template <int dim>
void FieldCurrentsHeating<dim>::setup_dofs () {
	set_active_fe_indices ();
	dof_handler.distribute_dofs (fe_collection);
	{
		constraints.clear ();
		// Actually not using mesh refinement (or hanging nodes)
		//DoFTools::make_hanging_node_constraints (dof_handler, constraints);

		/* Apply Dirichlet boundary conditions in assemble_system so that
		 * the boundaries could be changed depending on newton iteration number

		const FEValuesExtractors::Scalar potential_c(1);
		VectorTools::interpolate_boundary_values (dof_handler,
												  BoundaryId::bulk_bottom,
												  ZeroFunction<dim>(3),
												  constraints,
												  fe_collection.component_mask(potential_c));
		const FEValuesExtractors::Scalar temperature(2);
		VectorTools::interpolate_boundary_values (dof_handler,
												  BoundaryId::bulk_bottom,
												  ConstantFunction<dim>(300.0, 3),
												  constraints,
												  fe_collection.component_mask(temperature));
		*/
	}

	{
		// Set vacuum potential to be 0 at the copper surface (internal boundary)

		std::vector<types::global_dof_index> local_face_dof_indices (vacuum_fe.dofs_per_face);
		for (typename hp::DoFHandler<dim>::active_cell_iterator
				cell = dof_handler.begin_active();
				cell != dof_handler.end(); ++cell) {
			if (cell_is_in_vacuum_domain (cell)) {
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
					if (!cell->at_boundary(f) && cell_is_in_copper_domain(cell->neighbor(f))) {
						cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
						for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
							if (vacuum_fe.face_system_to_component_index(i).first < dim)
								constraints.add_line (local_face_dof_indices[i]);
					}
				}
			}
		}
	}

	constraints.close ();
	std::cout << "   Number of active cells: "
			  << triangulation.n_active_cells()
			  << std::endl
			  << "   Number of degrees of freedom: "
			  << dof_handler.n_dofs()
			  << std::endl;
	{
		// Create sparsity pattern so that vacuum potential couples to itself in cells
		// copper variables couple to each other in cells
		// and vacuum potential couples to copper potential on faces
		// NB: actually vac_pot doesn't need to couple to cop_pot and is commented out

		DynamicSparsityPattern dsp (dof_handler.n_dofs(),
									dof_handler.n_dofs());
		Table<2,DoFTools::Coupling> cell_coupling (fe_collection.n_components(),
												   fe_collection.n_components());
		Table<2,DoFTools::Coupling> face_coupling (fe_collection.n_components(),
												   fe_collection.n_components());

		// pot_vacuum = 0, pot_copper = 1, temperature = 2
		for (unsigned int c=0; c<fe_collection.n_components(); ++c)
			for (unsigned int d=0; d<fe_collection.n_components(); ++d) {
				if (c == d || (c == 1 && d == 2) || (c == 2 && d == 1))
					cell_coupling[c][d] = DoFTools::always;
				//if ((c == 0 && d == 1) || (c == 1 && d == 0))
				//	face_coupling[c][d] = DoFTools::always;
			}
		DoFTools::make_flux_sparsity_pattern (dof_handler, dsp,
											  cell_coupling, face_coupling);
		constraints.condense (dsp);
		sparsity_pattern.copy_from (dsp);
	}

	system_matrix.reinit (sparsity_pattern);
	solution.reinit (dof_handler.n_dofs());
	previous_solution.reinit (dof_handler.n_dofs());
	system_rhs.reinit (dof_handler.n_dofs());
}

template <int dim>
void FieldCurrentsHeating<dim>::assemble_system_newton (const bool first_iteration) {
	system_matrix=0;
	system_rhs=0;
	const QGauss<dim> vacuum_quadrature(field_degree+2);
	// Not a good idea to have different quadratures for currents and heating (probably?)
	const QGauss<dim> copper_quadrature(std::max(currents_degree, heating_degree)+2);
	hp::QCollection<dim>  q_collection;
	q_collection.push_back (vacuum_quadrature);
	q_collection.push_back (copper_quadrature);
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection,
								  update_values    |
								  update_quadrature_points  |
								  update_JxW_values |
								  update_gradients);
	const QGauss<dim-1> common_face_quadrature(std::max (field_degree+2,
															currents_degree+2));
	FEFaceValues<dim>    field_fe_face_values (vacuum_fe,
											  common_face_quadrature,
											  update_values |
											  update_JxW_values |
											  update_normal_vectors |
											  update_gradients);
	FEFaceValues<dim>    copper_fe_face_values (copper_fe,
												  common_face_quadrature,
												  update_values |
												  update_JxW_values |
												  update_normal_vectors |
												  update_gradients);

	const unsigned int        vacuum_dofs_per_cell  = vacuum_fe.dofs_per_cell;
	const unsigned int        copper_dofs_per_cell 	= copper_fe.dofs_per_cell;
	FullMatrix<double>        local_matrix;
	FullMatrix<double>        local_interface_matrix (copper_dofs_per_cell,
													  copper_dofs_per_cell);
	Vector<double>            local_rhs;
	Vector<double>            local_interface_rhs (copper_dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices;
	std::vector<types::global_dof_index> neighbor_dof_indices (vacuum_dofs_per_cell);

	// ---------------------------------------------------------------------------------------------
	// Declaring variables once at the start (for performance reasons, probably?)

	// The previous solution values in the cell quadrature points
	std::vector<Tensor<1,dim>> prev_sol_potential_gradients (copper_quadrature.size());
	std::vector<double> prev_sol_temperature_values (copper_quadrature.size());
	std::vector<Tensor<1,dim>> prev_sol_temperature_gradients (copper_quadrature.size());

	// The previous solution values in the face quadrature points
	std::vector<Tensor<1,dim>> prev_sol_face_potential_gradients (common_face_quadrature.size());
	std::vector<Tensor<1,dim>> prev_sol_face_temperature_gradients (common_face_quadrature.size());
	std::vector<double> prev_sol_face_temperature_values (common_face_quadrature.size());
	std::vector<Tensor<1,dim>> prev_sol_face_field (common_face_quadrature.size());

	// Shape function values and gradients (arrays for every cell DOF)
	std::vector< Tensor<1, dim> >	potential_c_phi_grad (copper_dofs_per_cell);
	std::vector<double> 			temperature_phi (copper_dofs_per_cell);
	std::vector< Tensor<1, dim> >	temperature_phi_grad (copper_dofs_per_cell);
	// ----------------------------------------------------------------------------------------------

	const FEValuesExtractors::Scalar     potential_v (0);
	const FEValuesExtractors::Scalar     potential_c (1);
	const FEValuesExtractors::Scalar     temperature (2);

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell) {
	    hp_fe_values.reinit (cell);
	    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
	    local_matrix.reinit (cell->get_fe().dofs_per_cell,
	                         cell->get_fe().dofs_per_cell);
	    local_rhs.reinit (cell->get_fe().dofs_per_cell);

	    // -----------------------------------------------------------------------------------------------
	    // Cell is in vacuum domains
	    // -------------------------------
	    if (cell_is_in_vacuum_domain (cell)) {
	    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	    	Assert (dofs_per_cell == vacuum_dofs_per_cell, ExcInternalError());

	    	// Integration (?) over the cell -> create matrix entries
	        for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q) {
	            for (unsigned int i=0; i<dofs_per_cell; ++i)
	              for (unsigned int j=0; j<dofs_per_cell; ++j)
	                local_matrix(i,j) += fe_values[potential_v].gradient(i, q)
	                					 * fe_values[potential_v].gradient(j, q)
	                                     * fe_values.JxW(q);
	        }
	        // Integration over cell's boundaries (if there is Neumann BC)
	        // Only for 1st Newton iteration
	        if (first_iteration) {
				for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
					// If face is on vacuum top boundary, apply Neumann BC
					if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == BoundaryId::vacuum_top) {
						field_fe_face_values.reinit(cell, f);

						for (unsigned int q_index = 0; q_index < field_fe_face_values.n_quadrature_points; ++q_index) {
							for (unsigned int i = 0; i < dofs_per_cell; ++i) {
								local_rhs(i) += (field_fe_face_values[potential_v].value(i, q_index)
										* applied_field * field_fe_face_values.JxW(q_index));
							}
						}
					}
				}
	        }
		// -----------------------------------------------------------------------------------------------
		// Cell is in copper domain
		// -------------------------------
	    } else {
	        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	        Assert (dofs_per_cell == copper_dofs_per_cell, ExcInternalError());

	        fe_values[potential_c].get_function_gradients(previous_solution, prev_sol_potential_gradients);
			fe_values[temperature].get_function_values(previous_solution, prev_sol_temperature_values);
			fe_values[temperature].get_function_gradients(previous_solution, prev_sol_temperature_gradients);

			// Integration (?) over the cell -> create matrix entries
	        for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q) {

	        	double prev_temp = prev_sol_temperature_values[q];
	        	const Tensor<1, dim> prev_pot_grad = prev_sol_potential_gradients[q];
	        	const Tensor<1, dim> prev_temp_grad = prev_sol_temperature_gradients[q];

				double sigma = pq.sigma(prev_temp);
				double dsigma = pq.dsigma(prev_temp);
				double kappa = pq.kappa(prev_temp);
				double dkappa = pq.dkappa(prev_temp);

	            for (unsigned int k=0; k<dofs_per_cell; ++k) {
	            	potential_c_phi_grad[k] = fe_values[potential_c].gradient (k, q);
	            	temperature_phi[k] = fe_values[temperature].value (k, q);
	            	temperature_phi_grad[k] = fe_values[temperature].gradient (k, q);
	            }
	            for (unsigned int i=0; i<dofs_per_cell; ++i) {
	            	for (unsigned int j=0; j<dofs_per_cell; ++j) {
	            		local_matrix(i, j)
							+= ( - (potential_c_phi_grad[i] * sigma * potential_c_phi_grad[j])
								 - (potential_c_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
								 + (temperature_phi[i] * 2 * sigma * prev_pot_grad * potential_c_phi_grad[j])
								 + (temperature_phi[i] * dsigma * prev_pot_grad * prev_pot_grad * temperature_phi[j])
								 - (temperature_phi_grad[i] * dkappa * prev_temp_grad * temperature_phi[j])
								 - (temperature_phi_grad[i] * kappa * temperature_phi_grad[j])
								 ) * fe_values.JxW(q);
	            	}
	            	local_rhs(i)
	            		+= (   (potential_c_phi_grad[i] * sigma * prev_pot_grad)
	            			 - (temperature_phi[i] * sigma * prev_pot_grad * prev_pot_grad)
							 + (temperature_phi_grad[i] * kappa * prev_temp_grad)
	            			) * fe_values.JxW(q);
	            }

	        }
	    }

	    // ---------------------------------------------------------------------------------------------

	    local_dof_indices.resize (cell->get_fe().dofs_per_cell);
	    cell->get_dof_indices (local_dof_indices);
	    constraints.distribute_local_to_global (local_matrix, local_rhs,
	                                            local_dof_indices,
	                                            system_matrix, system_rhs);

	    // ---------------------------------------------------------------------------------------------
	    // Take care of the boundary between vacuum and copper
	    if (cell_is_in_copper_domain(cell)) {
	    	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
	    		if (cell->at_boundary(f) == false) {
	    			if (cell_is_in_vacuum_domain(cell->neighbor(f))) {
	    				copper_fe_face_values.reinit(cell, f);
	    				field_fe_face_values.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));

	    				// "assemble_interface_term" starts
	    				Assert (field_fe_face_values.n_quadrature_points == copper_fe_face_values.n_quadrature_points,
	    				          ExcInternalError());
	    				const unsigned int n_face_quadrature_points = field_fe_face_values.n_quadrature_points;

	    				local_interface_matrix = 0;
	    				local_interface_rhs = 0;

	    				copper_fe_face_values[potential_c].get_function_gradients(previous_solution, prev_sol_face_potential_gradients);
	    				copper_fe_face_values[temperature].get_function_gradients(previous_solution, prev_sol_face_temperature_gradients);
	    				copper_fe_face_values[temperature].get_function_values(previous_solution, prev_sol_face_temperature_values);
	    				field_fe_face_values[potential_v].get_function_gradients(previous_solution, prev_sol_face_field);

	    				for (unsigned int q = 0; q < n_face_quadrature_points; ++q) {

	    					double prev_temp = prev_sol_face_temperature_values[q];
							double e_field = prev_sol_face_field[q].norm();
							double emission_current = pq.emission_current(e_field, prev_temp);

							// Nottingham heat flux in
							// (eV*A/nm^2) -> (eV*n*q_e/(s*nm^2)) -> (J*n/(s*nm^2)) -> (W/nm^2)
							double nottingham_flux = pq.nottingham_de(e_field, prev_temp)*emission_current;
							nottingham_flux = 0.0;

							const Tensor<1,dim> prev_face_pot_grad = prev_sol_face_potential_gradients[q];
							const Tensor<1,dim> prev_face_temp_grad = prev_sol_face_temperature_gradients[q];

							//double sigma = pq.sigma(prev_temp);
							//double kappa = pq.kappa(prev_temp);
							double dsigma = pq.dsigma(prev_temp);
							double dkappa = pq.dkappa(prev_temp);

							const Tensor<1,dim> face_normal_vector = copper_fe_face_values.normal_vector(q);

							for (unsigned int i = 0; i < copper_fe_face_values.dofs_per_cell; ++i) {
								const double potential_c_phi_i = copper_fe_face_values[potential_c].value(i, q);
								const double temperature_phi_i = copper_fe_face_values[temperature].value(i, q);

								local_interface_rhs(i) += ( - potential_c_phi_i * emission_current
															- temperature_phi_i * nottingham_flux
															) * copper_fe_face_values.JxW(q);

								/*
								local_interface_rhs(i) += ( - potential_c_phi_i * face_normal_vector * prev_face_pot_grad * sigma
															- temperature_phi_i * face_normal_vector * prev_face_temp_grad * kappa
															) * copper_fe_face_values.JxW(q);
								*/
								for (unsigned int j = 0; j < copper_fe_face_values.dofs_per_cell; ++j) {
									const double temperature_phi_j = copper_fe_face_values[temperature].value(j, q);

									local_interface_matrix(i, j)
										+= ( potential_c_phi_i * face_normal_vector * dsigma
												* prev_face_pot_grad * temperature_phi_j
											+ temperature_phi_i * face_normal_vector * dkappa
												* prev_face_temp_grad * temperature_phi_j
										   ) * copper_fe_face_values.JxW(q);

								}
							}

	    				}
	    				// "assemble_interface_term" ends
	    				cell->neighbor(f)->get_dof_indices (neighbor_dof_indices);
	    				constraints.distribute_local_to_global (local_interface_matrix, local_interface_rhs,
																local_dof_indices,
																system_matrix, system_rhs);

	    			}
	    		}
	    	}
	    }
	    // ---------------------------------------------------------------------------------------------

	}

	// ----------------------------------------------------------------------------------------------------------
	// Apply Dirichlet boundary conditions
	// Note: if it was done outside assembly, e.g. in setup_dofs using constraints, performance would be better

	// 0 potential at bulk bottom
	std::map<types::global_dof_index, double> potential_c_dirichlet;
	VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::bulk_bottom, ZeroFunction<dim>(3),
			potential_c_dirichlet, fe_collection.component_mask(potential_c));
	MatrixTools::apply_boundary_values(potential_c_dirichlet, system_matrix,
			solution, system_rhs);

	// 300K at bulk bottom if initial step, 0 otherwise
	std::map<types::global_dof_index, double> temperature_dirichlet;
	if (first_iteration) {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::bulk_bottom,
				ConstantFunction<dim>(300.0, 3), temperature_dirichlet, fe_collection.component_mask(temperature));
	} else {
		VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::bulk_bottom,
				ZeroFunction<dim>(3), temperature_dirichlet, fe_collection.component_mask(temperature));
	}
	MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix,
			solution, system_rhs);
	// ----------------------------------------------------------------------------------------------------------
}


template <int dim>
void FieldCurrentsHeating<dim>::solve () {
	SparseDirectUMFPACK direct_solver;
	direct_solver.initialize (system_matrix);
	direct_solver.vmult (solution, system_rhs);
	constraints.distribute (solution);
}

// ----------------------------------------------------------------------------------------
// Class for outputting the resulting field distribution (calculated from potential distr.)
template <int dim>
class EFieldPostProcessor : public DataPostprocessorVector<dim> {
public:
	EFieldPostProcessor() : DataPostprocessorVector<dim>("electric_field", update_gradients) {}
	void
	compute_derived_quantities_vector (	const std::vector< Vector< double > > &  					/*uh*/,
											const std::vector< std::vector< Tensor< 1, dim > > > &  duh,
											const std::vector< std::vector< Tensor< 2, dim > > > &  /*dduh*/,
											const std::vector< Point< dim > > &  					/*normals*/,
											const std::vector< Point< dim > > &  					/*evaluation_points*/,
											std::vector< Vector< double > > &  						computed_quantities
										  ) const {
		for (unsigned int i=0; i<computed_quantities.size(); i++) {
			for (unsigned int d=0; d<dim; ++d)
				computed_quantities[i](d) = duh[i][0][d];
		}
	}
};
// Class for outputting the resulting current distribution
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
			double t = uh[i][2]; // temperature
			double sigma = pq.sigma(t);
			for (unsigned int d=0; d<dim; ++d) {
				double e_field = -duh[i][1][d]; // 0 - pot_v, 1 - pot_c, 2 - temp
				computed_quantities[i](d) = sigma*e_field;
			}
		}
	}
};
// ----------------------------------------------------------------------------

template <int dim>
void FieldCurrentsHeating<dim>::output_results (const unsigned int iteration) const {
	std::string filename = "solution-fch-" + Utilities::int_to_string(iteration) + ".vtk";

	std::vector<std::string> solution_names {"potential_v", "potential_c", "temperature"};

	EFieldPostProcessor<dim> efield_post_processor; // needs to be before data_out
	CurrentPostProcessor<dim> current_post_processor(pq); // needs to be before data_out

	DataOut<dim,hp::DoFHandler<dim> > data_out;
	data_out.attach_dof_handler (dof_handler);

	data_out.add_data_vector(solution, solution_names);
	data_out.add_data_vector(solution, efield_post_processor);
	data_out.add_data_vector(solution, current_post_processor);
	data_out.build_patches();

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}


template <int dim>
void FieldCurrentsHeating<dim>::run () {
	Timer timer;
	timer.start ();
	setup_dofs();
	timer.stop ();
	deallog << "setup_dofs: " << timer () << "s" << std::endl;


	for (unsigned int it = 0; it < 10; it++) {
		std::cout << "    Newton iteration " << it << std::endl;

		system_matrix.reinit (sparsity_pattern);
		system_rhs.reinit(dof_handler.n_dofs());


		timer.restart ();
		assemble_system_newton(it == 0);
		timer.stop ();

		deallog << "system assembled in "
				<< timer ()
				<< "s"
				<< std::endl;

		timer.restart();
		solve();
		solution.add(1.0, previous_solution);
		timer.stop();

		deallog << "solver finished in "
				<< timer ()
				<< "s"
				<< std::endl;

		output_results(it);

		previous_solution = solution;
	}
}

template class FieldCurrentsHeating<2> ;
template class FieldCurrentsHeating<3> ;

} // namespace field_currents_heating
