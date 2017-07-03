/*
 * laplace.h
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_LAPLACE_H_
#define INCLUDE_LAPLACE_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>
#include <iostream>

#include "currents_and_heating.h" // for friend class declaration
#include "currents_and_heating_stationary.h" // for friend class declaration
#include "mesh_preparer.h"

namespace fch {

// forward declaration for friend class
template<int dim> class CurrentsAndHeating;
template<int dim> class CurrentsAndHeatingStationary;

using namespace dealii;

template<int dim>
class Laplace {
public:
    Laplace();

    /** runs the calculation */
    void run();

    /** getter for the mesh */
    Triangulation<dim>* get_triangulation();
    /** getter for dof_handler */
    DoFHandler<dim>* get_dof_handler();

    /** Sets the applied electric field in GV/m (V/nm) boundary condition */
    void set_applied_efield(const double applied_field_);

    /**
     * Imports mesh from file and optionally outputs it to a .vtk file
     * Additionally sets the boundary indicators corresponding to vacuum
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
    void output_mesh(const std::string file_name = "vacuum_mesh.vtk");

    /** get the electric field at the specified point */
    double probe_efield(const Point<dim> &p) const;

    /**
     * method to obtain the electric potential values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return potential values in the specified nodes
     */
    std::vector<double> get_potential(const std::vector<int> &cell_indexes,
            const std::vector<int> &vert_indexes);

    /**
     * method to obtain the electric field values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return electric field vectors in the specified nodes
     */
    std::vector<Tensor<1, dim>> get_efield(const std::vector<int> &cell_indexes,
            const std::vector<int> &vert_indexes);

    /** sets up the system: calculates the sparse matrix dimensions, etc */
    void setup_system();
    /** assembles the matrix equation */
    void assemble_system();

    /** solves the matrix equation using conjugate gradients
     * @param max_iter maximum number of iterations allowed
     * @param tol tolerance
     * @param pc_ssor flag to use SSOR preconditioner
     * @param ssor_param parameter to SSOR preconditioner, 1.2 is known to work well with laplace
     */
    void solve(int max_iter = 2000, double tol = 1e-9, bool pc_ssor = true,
            double ssor_param = 1.2);

    /** outputs the results to a specified file */
    void output_results(
            const std::string filename = "field_solution.vtk") const;
    /** equivalent to output_results */
    void write(const std::string filename = "field_solution.vtk") const {
        output_results(filename);
    };

    /** string stream prints the statistics about the system */
    friend std::ostream& operator <<(std::ostream &os, const Laplace<dim>& d) {
        os << "#elems=" << d.triangulation.n_active_cells() << ",\t#faces="
                << d.triangulation.n_active_faces() << ",\t#edges="
                << d.triangulation.n_active_lines() << ",\t#nodes="
                << d.triangulation.n_used_vertices() << ",\t#dofs="
                << d.dof_handler.n_dofs();
        return os;
    }

private:
    static constexpr unsigned int shape_degree = 1;
    static constexpr unsigned int quadrature_degree = shape_degree + 1;

    static constexpr double applied_efield_default = 2.0;
    double applied_efield;

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    friend class CurrentsAndHeating<dim> ;
    friend class CurrentsAndHeatingStationary<dim> ;
};

} // namespace fch

#endif /* INCLUDE_LAPLACE_H_ */
