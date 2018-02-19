/*
 * laplace.h -> Laplace.h
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan
 */

#ifndef LAPLACE_H_
#define LAPLACE_H_

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

#include "CurrentsAndHeating.h" // for friend class declaration
#include "CurrentsAndHeatingStationary.h" // for friend class declaration
#include "MeshPreparer.h"
#include "ParticleSpecies.h"


#include "DealSolver.h"

namespace fch {

// forward declaration for friend class
template<int dim> class CurrentsAndHeating;
template<int dim> class CurrentsAndHeatingStationary;
template<int dim> class CurrentHeatStatSolver;


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

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p) const;

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p, int cell_index) const;

    /** get the potential value at a specified point using dealii (slow)  */
    double probe_potential(const Point<dim> &p) const;

    /** get the potential value at a specified point using dealii with known cell id for the point  */
    double probe_potential(const Point<dim> &p, const int cell_index) const;

    /** Probes the field at point p that belongs in cell with cell_index. Fast, when cell_index is correct */
    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index) const;

    /** saves the system matrix (useful before BCs have been applied) */
    void save_system(){
        system_matrix_save.copy_from(system_matrix);
//        system_rhs_save = system_rhs;
    }

    /** Restores the system matrix to the saved one */
    void restore_system(){
        system_matrix.copy_from(system_matrix_save);
//        system_rhs = system_rhs_save;
    }

    /** Values of the shape functions at point p with respect to the nodes of
     * cell with cell_index p has to belong in cell_index!!
     */
    std::vector<double> shape_funs(const Point<dim> &p, int cell_index) const;
    /** Helper function for the above */
    std::vector<double> shape_funs(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /**
     * method to obtain the electric potential values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return potential values in the specified nodes
     */
    std::vector<double> get_potential(const std::vector<int> &cell_indexes,
            const std::vector<int> &vert_indexes);


    double get_cell_vol(int cellid);

    /**
     * method to obtain the electric field values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return electric field vectors in the specified nodes
     */
    std::vector<Tensor<1, dim>> get_efield(const std::vector<int> &cell_indexes,
            const std::vector<int> &vert_indexes) const;

    /** @brief set up dynamic sparsity pattern
     *  a) define optimal structure for sparse matrix representation,
     *  b) allocate memory for sparse matrix and solution and right-hand-side (rhs) vector
     */
    void setup_system(bool first_time);

    /** @brief assemble the matrix equation
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the Laplace equation weak formulation and to the boundary conditions.
     */
    void assemble_system(double applied_field){
        assemble_system_lhs();
        assemble_system_neuman(BoundaryId::vacuum_top, applied_field);
        assemble_system_dirichlet(BoundaryId::copper_surface, 0.0);
        assemble_system_finalize();
    }

    /** @brief Reset the system and assemble the LHS matrix
     * Calculate sparse matrix elements
     * according to the Laplace equation weak formulation
     * This should be the first function call to setup the equations (after setup_system() ).
     */
    void assemble_system_lhs();

    /** @brief Initialization of the RHS of the matrix equation
     * Set the right-hand-side vector for Neuman boundary conditions on the given BoundaryId.
     */
    void assemble_system_neuman(BoundaryId bid, double E0);

    /** @brief Assemble the RHS of the matrix equation
     * Add to the right-hand-side vector for point charges, as used in PIC.
     */
    void assemble_system_pointcharge(femocs::ParticleSpecies &parts);

    /** @brief Give value potential to all DOFs with BoundaryId bid
     */
    void assemble_system_dirichlet(BoundaryId bid, double potential);

    /** @brief Apply all dirichlet boundary conditions to the system.
     * This should be the last function call to setup the equations, before calling solve()
     */
    void assemble_system_finalize() {
        MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
    }

    /** solves the matrix equation using conjugate gradient method
     * @param max_iter maximum number of iterations allowed
     * @param tol tolerance of the solution
     * @param pc_ssor flag to use SSOR preconditioner
     * @param ssor_param   parameter to SSOR preconditioner. 1.2 is known to work well with laplace.
     *                     its fine tuning optimises calculation time
     */
    int solve(int max_iter = 2000, double tol = 1e-9, bool pc_ssor = true,
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

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;     ///< structure for sparse matrix representation
    SparseMatrix<double> system_matrix;   ///< system matrix of matrix equation
    SparseMatrix<double> system_matrix_save;   ///< system matrix of matrix equation (save before Dirichlet BCs remove dofs)

    Vector<double> solution;              ///< resulting electric potential in the mesh nodes
    Vector<double> system_rhs;            ///< right-hand-side of the matrix equation
    Vector<double> system_rhs_save; ///< system RHS of matrix equation (save before Dirichlet BCs remove dofs)

    std::map<types::global_dof_index, double> boundary_values; // Map of dirichlet boundary conditions

    double probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    double probe_efield_norm(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    friend class CurrentsAndHeating<dim> ;
    friend class CurrentsAndHeatingStationary<dim> ;
};



using namespace std;

/** @brief Class to solve Laplace equation in 2D or 3D
 * It is inspired by the step-3 of Deal.II tutorial
 * https://www.dealii.org/8.5.0/doxygen/deal.II/step_3.html
 */
template<int dim>
class PoissonSolver : public DealSolver<dim> {
public:
    PoissonSolver();

    /** Run Conjugate-Gradient solver to solve matrix equation */
    unsigned int solve() { return this->solve_cg(2000, 1e-9, 1.2); }
    
    /** Runs the calculation: setup and assemble system, solve PoissonSolver equation, output the results*/
    void run();

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p) const;

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p, int cell_index) const;

    /** get the potential value at a specified point using dealii (slow)  */
    double probe_potential(const Point<dim> &p) const;

    /** get the potential value at a specified point using dealii with known cell id for the point  */
    double probe_potential(const Point<dim> &p, const int cell_index) const;

    /** Probes the field at point p that belongs in cell with cell_index. Fast, when cell_index is correct */
    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index) const;

    /**
     * method to obtain the electric potential values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return potential values in the specified nodes
     */
    vector<double> get_potential(const vector<int> &cell_indexes,
            const vector<int> &vert_indexes);

    /**
     * method to obtain the electric field values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return electric field vectors in the specified nodes
     */
    vector<Tensor<1, dim>> get_efield(const vector<int> &cell_indexes,
            const vector<int> &vert_indexes) const;

    /** @brief set up dynamic sparsity pattern
     *  a) define optimal structure for sparse matrix representation,
     *  b) allocate memory for sparse matrix and solution and right-hand-side (rhs) vector
     */
    void setup_system(bool first_time);

    /** @brief assemble the matrix equation
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the PoissonSolver equation weak formulation and to the boundary conditions.
     */
    void assemble_system(double applied_field);

    /** Store long range electric field value */
    void set_applied_field(const double applied_field_) { applied_field = applied_field_; }

    /** @brief Reset the system and assemble the LHS matrix
     * Calculate sparse matrix elements
     * according to the Laplace equation weak formulation
     * This should be the first function call to setup the equations (after setup_system() ).
     */
    void assemble_lhs();

    /** @brief Assemble the RHS of the matrix equation
     * Add to the right-hand-side vector for point charges, as used in PIC.
     */
    void assemble_system_pointcharge(femocs::ParticleSpecies &parts);

    /** @brief Give value potential to all DOFs with BoundaryId bid
     */
    void assemble_system_dirichlet(BoundaryId bid, double potential);

    /** @brief Apply all dirichlet boundary conditions to the system.
     * This should be the last function call to setup the equations, before calling solve()
     */
    void assemble_system_finalize();

    /** Outputs the results (electric potential and field) to a specified file in vtk format */
    void output_results(const string &filename) const;

private:
    map<types::global_dof_index, double> boundary_values; // Map of dirichlet boundary conditions
    double applied_field; ///< long-range electric field

    double probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;
    
    double probe_efield_norm(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Mark different regions in mesh */
    void mark_boundary();

    /** Return the boundary condition value at the centroid of face */
    double get_face_bc(const unsigned int face) const;

    friend class CurrentsAndHeating<dim> ;
    friend class CurrentsAndHeatingStationary<dim> ;
    friend class CurrentHeatStatSolver<dim> ;
};

} // namespace fch

#endif /* LAPLACE_H_ */
