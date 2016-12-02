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

#include "mesh_preparer.h" // for BoundaryId-s.. probably should think of a better place for them
#include "physical_quantities.h"
#include "laplace.h"


namespace fch {

	// forward declaration for Laplace to exist when declaring CurrentsAndHeating
	template<int dim> class Laplace;

	using namespace dealii;

	template<int dim>
	class CurrentsAndHeating {
	public:

		/**
		 * Initializes the object
		 * NB: pq and laplace will be set as NULL, so they must be set separately
		 */
		CurrentsAndHeating();

		/**
		 * Initializes the object without initial condition interpolation
		 * @param pq_ object from physical data (emission currents, resistance, etc) is obtained
		 * @param laplace_ electric field solver for the corresponding vacuum domain
		 */
		CurrentsAndHeating(PhysicalQuantities *pq_, Laplace<dim> *laplace_);

		/**
		 * Initializes the object with initial condition interpolation from another CurrentsAndHeating object
		 * @param pq_ object from physical data (emission currents, resistance, etc) is obtained
		 * @param laplace_ electric field solver for the corresponding vacuum domain
		 * @param ch_previous_iteration_ object from where the initial condition is interpolated;
		 * 								 it must contain a solution and its mesh can be different than the current mesh
		 */
		CurrentsAndHeating(PhysicalQuantities *pq_, Laplace<dim> *laplace_,
						   CurrentsAndHeating *ch_previous_iteration_);

		/**
		 * Reinitializes current object without initial condition interpolation
		 * The mesh must be imported again (corresponding to the new laplace object)
		 */
		void reinitialize(Laplace<dim> *laplace_);

		/**
		 * Reinitializes current object with initial condition interpolation
		 * The mesh must be imported again (corresponding to the new laplace object)
		 */
		void reinitialize(Laplace<dim> *laplace_, CurrentsAndHeating *ch_previous_iteration_);

		/** Sets the physical quantities object */
		void set_physical_quantities(PhysicalQuantities *pq_);

		/** Sets the ambient temperature boundary condition */
		void set_ambient_temperature(const double ambient_temperature_);

		/** runs the calculation with hardcoded parameters (mainly for testing) */
		void run();

		/**
		 * Runs calculation with specific parameters
		 * @param temperature_tolerance the tolerance of temperature when newton iterations are stopped
		 * @param max_newton_iter maximum number of newton iterations
		 * @param file_output should the solution be output to a .vtk file
		 * @param out_fname if file_output is set to true, then the newton iterations are saved to files <out_fname>-<N#>.vtk,
		 * 					where N# is the number of the newton iteration
		 * @param print boolean if calculation info should be output to cout
		 */
		double run_specific(double temperature_tolerance=1.0, int max_newton_iter=10,
						  bool file_output=true, std::string out_fname="sol", bool print=true);

		/** getter for the mesh */
		Triangulation<dim>* get_triangulation();

		/** getter for dof_handler */
		DoFHandler<dim>* get_dof_handler();

		/** getter for the solution */
		Vector<double>* get_solution();

		/**
		 * Imports mesh from file and optionally outputs it to a .vtk file
		 * Additionally sets the boundary indicators corresponding to copper
		 * @param file_name file from the mesh is imported
		 * @param out_name if empty, won't save the file, otherwise saves the mesh to .vtk file
		 */
		void import_mesh_from_file(const std::string file_name, const std::string out_name = "");

		/**
		 * imports mesh directly from vertex and cell data and sets the boundary indicators
		 * @return true if success, otherwise false
		 */
		bool import_mesh_directly(std::vector<Point<dim> > vertices, std::vector<CellData<dim> > cells);

		/** outputs the mesh to .vtk file */
		void output_mesh(const std::string file_name = "copper_mesh.vtk");

		/**
		 * method to obtain the temperature values in selected nodes
		 * @param cell_indexes global cell indexes, where the corresponding nodes are situated
		 * @param vert_indexes the vertex indexes of the nodes inside the cell
		 * @return temperature values in the specified nodes
		 */
		std::vector<double> get_temperature(const std::vector<int> &cell_indexes,
										  	const std::vector<int> &vert_indexes);
		/**
		 * method to obtain the current density values in selected nodes
		 * @param cell_indexes global cell indexes, where the corresponding nodes are situated
		 * @param vert_indexes the vertex indexes of the nodes inside the cell
		 * @return current density vectors in the specified nodes
		 */
		std::vector<Tensor<1, dim>> get_current(const std::vector<int> &cell_indexes,
									  	  	  	const std::vector<int> &vert_indexes);

	    /** string stream prints the statistics about the system */
	    friend std::ostream& operator <<(std::ostream &os, const CurrentsAndHeating<dim>& d) {
	        os << "#elems=" << d.triangulation.n_active_cells()
	                << ",\t#faces=" << d.triangulation.n_active_faces()
	                << ",\t#edges=" << d.triangulation.n_active_lines()
	                << ",\t#nodes=" << d.triangulation.n_used_vertices()
	                << ",\t#dofs=" << d.dof_handler.n_dofs();
	        return os;
	    }

	private:
		void setup_system();
		void assemble_system_newton();
		void solve();
		void output_results(const unsigned int iteration, const std::string fname = "sol") const;

		void setup_mapping();
		void set_initial_condition();
		void set_initial_condition_slow();

		static constexpr unsigned int currents_degree = 1;
		static constexpr unsigned int heating_degree  = 1;

		static constexpr double ambient_temperature_default = 300.0;
		double ambient_temperature;

		FESystem<dim> fe;

		Triangulation<dim> triangulation;
		DoFHandler<dim> dof_handler;

		SparsityPattern sparsity_pattern;
		SparseMatrix<double> system_matrix;

		Vector<double> present_solution;
		Vector<double> newton_update;
		Vector<double> system_rhs;

		PhysicalQuantities *pq;
		Laplace<dim> *laplace;

		/** Bijective mapping of interface faces
		 * (copper_cell_index, copper_cell_face) <-> (vacuum_cell_index, vacuum_cell_face)
		 */
		std::map< std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > interface_map;

		/** Previous iteration mesh and solution for setting the initial condition */
		CurrentsAndHeating* previous_iteration;
		bool interp_initial_conditions;
	};

} // end fch namespace


#endif /* INCLUDE_CURRENTS_AND_HEATING_H_ */
