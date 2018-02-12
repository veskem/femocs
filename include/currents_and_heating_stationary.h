/*
 * currents_and_heating.h
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_CURRENTS_AND_HEATING_STATIONARY_H_
#define INCLUDE_CURRENTS_AND_HEATING_STATIONARY_H_

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

/** @brief Class to solve stationary heat and continuity equations in 2D or 3D to obtain temperature and current density distribution in material.
 * The information about solving non-linear equations can be found in section 31 of lecture series on
 * http://www.math.colostate.edu/~bangerth/videos.html
 *
 * The Newton's method is covered in several Deal.II tutorials, however they cover other equations than needed here.
 * https://www.dealii.org/8.5.0/doxygen/deal.II/Tutorial.html
 */
template<int dim>
class CurrentsAndHeatingStationary {
public:

    /**
     * Initializes the object
     * NB: pq and laplace will be set as NULL, so they must be set separately
     */
    CurrentsAndHeatingStationary();

    /**
     * Initializes the object without initial condition interpolation
     * @param pq_ object from physical data (emission currents, resistance, etc) is obtained
     * @param laplace_ electric field solver for the corresponding vacuum domain
     */
    CurrentsAndHeatingStationary(PhysicalQuantities *pq_, Laplace<dim> *laplace_);

    /**
     * Initializes the object with initial condition interpolation from another CurrentsAndHeating object
     * @param pq_ object from physical data (emission currents, resistance, etc) is obtained
     * @param laplace_ electric field solver for the corresponding vacuum domain
     * @param ch_previous_iteration_ object from where the initial condition is interpolated;
     * 								 it must contain a solution and its mesh can be different than the current mesh
     */
    CurrentsAndHeatingStationary(PhysicalQuantities *pq_, Laplace<dim> *laplace_,
            CurrentsAndHeatingStationary *ch_previous_iteration_);

    /**
     * Reinitializes current object without initial condition interpolation
     * The mesh must be imported again (corresponding to the new laplace object)
     */
    void reinitialize(Laplace<dim> *laplace_);

    /**
     * Reinitializes current object with initial condition interpolation
     * The mesh must be imported again (corresponding to the new laplace object)
     */
    void reinitialize(Laplace<dim> *laplace_,
            CurrentsAndHeatingStationary *ch_previous_iteration_);

    /** Sets up degrees of freedom and the sparsity pattern */
    void setup_system();

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
     * @param sor_alpha successive over-relaxation coefficient
     * @param ic_interp_treshold peak temperature value of the previous iteration, which determines if interpolation is done
     * @param skip_field_mapping skip the (cell face) <-> (field) mapping on the surface; the field BC must be set by other means
     * @return final temperature error
     */
    double run_specific(double temperature_tolerance = 1.0,
            int max_newton_iter = 10, bool file_output = true,
            std::string out_fname = "sol", bool print = true,
            double sor_alpha = 1.0, double ic_interp_treshold = 400,
            bool skip_field_mapping = false);

    /** Provide triangulation object to get access to the mesh data */
    Triangulation<dim>* get_triangulation();

    /** Provide dof_handler object to get access to the mesh data */
    DoFHandler<dim>* get_dof_handler();

    /** Return the solution vector with potential and temperature values */
    Vector<double>* get_solution();

    /**
     * Imports mesh from file and sets the boundary indicators corresponding to copper
     * @param file_name file from the mesh is imported
     */
    void import_mesh_from_file(const std::string file_name);

    /**
     * imports mesh directly from vertex and cell data and sets the boundary indicators
     * @return true if success, otherwise false
     */
    bool import_mesh_directly(std::vector<Point<dim> > vertices,
            std::vector<CellData<dim> > cells);

    /** outputs the mesh to .vtk file */
    void output_mesh(const std::string file_name = "copper_mesh.vtk");

    /**
     * Outputs the calculated solution to file
     * @param file_name name of the destination file
     * @param iteration if positive, will be appendend to the name before extension
     */
    void output_results(const std::string file_name = "sol.vtk",
            const int iteration = -1) const;

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
    std::vector<Tensor<1, dim>> get_current(
            const std::vector<int> &cell_indexes,
            const std::vector<int> &vert_indexes);

    /** string stream prints the statistics about the system */
    friend std::ostream& operator <<(std::ostream &os,
            const CurrentsAndHeatingStationary<dim>& d) {
        os << "#elems=" << d.triangulation.n_active_cells() << ",\t#faces="
                << d.triangulation.n_active_faces() << ",\t#edges="
                << d.triangulation.n_active_lines() << ",\t#nodes="
                << d.triangulation.n_used_vertices() << ",\t#dofs="
                << d.dof_handler.n_dofs();
        return os;
    }

    /** export the centroids of surface faces */
    void get_surface_nodes(std::vector<Point<dim>>& nodes);

    /** read the electric field norm on the centroids of surface faces */
    void set_electric_field_bc(const std::vector<double>& elfields);

private:
    void assemble_system_newton();
    void solve();

    bool setup_mapping();
    /**
     * @param smoothing replaces top given % by their average + stdev (if negative, will ignore)
     */
    bool setup_mapping_field(double smoothing = 0.01);

    void set_initial_condition();
    void set_initial_condition_slow();

    static constexpr unsigned int currents_degree = 1; ///< degree of the shape functions in current density solver
    static constexpr unsigned int heating_degree = 1;  ///< degree of the shape functions in temperature solver

    static constexpr double ambient_temperature_default = 300.0;  ///< default temperature applied on bottom of the material
    double ambient_temperature;

    static constexpr double temperature_stopping_condition = 1800.0;  ///< temperature after which the calculation is forced to stop

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

    /** Mapping of copper interface face to vacuum side
     * (copper_cell_index, copper_cell_face) <-> (vacuum_cell_index, vacuum_cell_face)
     */
    std::map<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > interface_map;
    /** Mapping of copper interface face to vacuum side e field norm
     * (copper_cell_index, copper_cell_face) <-> (electric field norm)
     */
    std::map<std::pair<unsigned, unsigned>, double> interface_map_field;

    /** Previous iteration mesh and solution for setting the initial condition */
    CurrentsAndHeatingStationary* previous_iteration;
    bool interp_initial_conditions;
};

} // end fch namespace

#endif /* INCLUDE_CURRENTS_AND_HEATING_STATIONARY_H_ */
