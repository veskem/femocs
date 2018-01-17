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

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>


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

/** @brief Class to solve Laplace equation in 2D or 3D
 * It is inspired by the step-3 of Deal.II tutorial
 * https://www.dealii.org/8.5.0/doxygen/deal.II/step_3.html
 */
template<int dim>
class Laplace {
public:
    Laplace();

    /** Runs the calculation: setup and assemble system, solve Laplace equation, output the results*/
    void run();

    /** getter for the mesh */
    Triangulation<dim>* get_triangulation();
    /** getter for dof_handler */
    DoFHandler<dim>* get_dof_handler();

    /** Sets the applied electric field in GV/m (V/nm) boundary condition */
    void set_applied_efield(const double applied_field_);

    /**
     * Imports mesh from file and sets the vacuum boundary indicators
     * @param file_name name of the mesh file
     */
    void import_mesh_from_file(const std::string file_name);

    /**
     * imports mesh directly from vertex and cell data and sets the boundary indicators
     * @return true if success, otherwise false
     */
    bool import_mesh_directly(std::vector<Point<dim> > vertices, std::vector<CellData<dim> > cells);

    /** outputs the mesh to file in .vtk format */
    void output_mesh(const std::string file_name = "vacuum_mesh.vtk");

    /** get the electric field at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    Tensor<1,dim> probe_efield(const Point<dim> &p) const;

    
    /** get the potential value at a specified point using dealii (slow)
     */
    double probe_potential(const Point<dim> &p) const;

    double probe_value(const Point<dim> &p, const int cell_index) const;

    double probe_value(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    std::vector<double> shape_funs(const Point<dim> &p, int cell_index) const;


    std::vector<double> shape_funs(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    void test_probe();
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

    /** @brief set up dynamic sparsity pattern
     *  a) define optimal structure for sparse matrix representation,
     *  b) allocate memory for sparse matrix and solution and right-hand-side (rhs) vector
     */
    void setup_system();

    /** @brief assemble the matrix equation
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the Laplace equation weak formulation and to the boundary conditions.
     */
    void assemble_system();

    /** solves the matrix equation using conjugate gradient method
     * @param max_iter maximum number of iterations allowed
     * @param tol tolerance of the solution
     * @param pc_ssor flag to use SSOR preconditioner
     * @param ssor_param   parameter to SSOR preconditioner. 1.2 is known to work well with laplace.
     *                     its fine tuning optimises calculation time
     */
    void solve(int max_iter = 2000, double tol = 1e-9, bool pc_ssor = true,
            double ssor_param = 1.2);

    /** Outputs the results (electric potential and field) to a specified file in vtk format */
    void output_results(const std::string filename = "field_solution.vtk") const;

    /** Outputs the results (electric potential and field) to a specified file in vtk format.
     * Identical to output_results */
    void write(const std::string filename = "field_solution.vtk") const {
        output_results(filename);
    };

    /** Print the statistics about the mesh and # degrees of freedom */
    friend std::ostream& operator <<(std::ostream &os, const Laplace<dim>& d) {
        os << "#elems=" << d.triangulation.n_active_cells() << ",\t#faces="
                << d.triangulation.n_active_faces() << ",\t#edges="
                << d.triangulation.n_active_lines() << ",\t#nodes="
                << d.triangulation.n_used_vertices() << ",\t#dofs="
                << d.dof_handler.n_dofs();
        return os;
    }

private:
    static constexpr unsigned int shape_degree = 1;   ///< degree of the shape functions (linear, quadratic etc elements)
    static constexpr unsigned int quadrature_degree = shape_degree + 1;  ///< degree of the Gaussian numerical integration

    static constexpr double applied_efield_default = 2.0;
    double applied_efield;

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;     ///< structure for sparse matrix representation
    SparseMatrix<double> system_matrix;   ///< system matrix of matrix equation

    Vector<double> solution;              ///< resulting electric potential in the mesh nodes
    Vector<double> system_rhs;            ///< right-hand-side of the matrix equation

    friend class CurrentsAndHeating<dim> ;
    friend class CurrentsAndHeatingStationary<dim> ;
};

} // namespace fch

#endif /* INCLUDE_LAPLACE_H_ */
