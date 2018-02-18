/*
 * DealSolver.h
 *
 *  Created on: 12.2.2018
 *      Author: veske
 */

#ifndef DEALSOLVER_H_
#define DEALSOLVER_H_

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/hp/fe_values.h>

#include "MeshPreparer.h"


using namespace dealii;
using namespace std;
namespace fch {

template<int dim> class PoissonSolver;
template<int dim> class CurrentHeatSolver;
template<int dim> class CurrentSolver;
template<int dim> class HeatSolver;

/** @brief General class to implement FEM solver in Deal.II
 */
template<int dim>
class DealSolver {
public:

    DealSolver();
    DealSolver(Triangulation<dim> *tria);
    virtual ~DealSolver() {}
    
    /** Solve the matrix equation using conjugate gradient method
     * @param n_steps     maximum number of iterations allowed
     * @param tolerance   tolerance of the solution
     * @param ssor_param  parameter to SSOR preconditioner. Its fine tuning optimises calculation time
     */
    unsigned int solve_cg(const int n_steps, const double tolerance, const double ssor_param);

    /** Provide triangulation object to get access to the mesh data */
    Triangulation<dim>* get_triangulation() { return &triangulation; }

    /** Provide dof_handler object to get access to the mesh data */
    DoFHandler<dim>* get_dof_handler() { return &dof_handler; }

    /** Return the solution vector */
    Vector<double>* get_solution() { return &solution; };

    /** Values of the shape functions at point p with respect to the nodes of
     * cell with cell_index p has to belong in cell_index!!
     */
    vector<double> shape_funs(const Point<dim> &p, int cell_index) const;

    /** Return the volume/area of i-th cell */
    double get_cell_vol(const int i) const;

    /** export the centroids of surface faces */
    void get_surface_nodes(vector<Point<dim>>& nodes) const;

    /**
     * Import mesh from file and set the boundary indicators corresponding to copper
     * @param file_name   file from the mesh is imported
     */
    void import_mesh_from_file(const string &file_name);

    /**
     * imports mesh directly from vertex and cell data and sets the boundary indicators
     * @return true if success, otherwise false
     */
    bool import_mesh_directly(vector<Point<dim>> vertices, vector<CellData<dim>> cells);

    /** outputs the mesh to .vtk file */
    void output_mesh(const string &file_name);

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

    /** Run the full calculation */
    void run() {}

    /** Get the solution value at the specified point. NB: Slow! */
    double probe_solution(const Point<dim> &p) const;

    /** Calculate max value from solution vector */
    double max_solution() const;

    /** Write the simulation results to file */
    virtual void output_results(const string &filename) const {};

    void write(const string &filename) const {
        output_results(filename);
    };

    /** Print the statistics about the mesh and # degrees of freedom */
    friend ostream& operator <<(ostream &os, const DealSolver<dim>& d) {
        os << "#elems=" << d.triangulation.n_active_cells()
                << ", #faces=" << d.triangulation.n_active_faces()
                << ", #edges=" << d.triangulation.n_active_lines()
                << ", #nodes=" << d.triangulation.n_used_vertices()
                << ", #dofs=" << d.dof_handler.n_dofs();
        return os;
    }

protected:
    static constexpr unsigned int shape_degree = 1;   ///< degree of the shape functions (linear, quadratic etc elements)
    static constexpr unsigned int quadrature_degree = shape_degree + 1;  ///< degree of the Gaussian numerical integration

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;     ///< structure for sparse matrix representation
    SparseMatrix<double> system_matrix;   ///< system matrix of matrix equation
    SparseMatrix<double> system_matrix_save;   ///< system matrix of matrix equation (save before Dirichlet BCs remove dofs)
    Vector<double> system_rhs;            ///< right-hand-side of the matrix equation
    Vector<double> system_rhs_save;       ///< system RHS of matrix equation (save before Dirichlet BCs remove dofs)

    Vector<double> solution;              ///< resulting solution in the mesh nodes
    Vector<double> solution_old;         ///< resulting solution in the mesh nodes

    /** Mark different regions of the mesh */
    virtual void mark_boundary() {}

    /** Helper function for the public shape_funs */
    vector<double> shape_funs(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    friend class PoissonSolver<dim>;
    friend class CurrentHeatSolver<dim>;
    friend class CurrentSolver<dim>;
    friend class HeatSolver<dim>;
};

} /* namespace femocs */

#endif /* DEALSOLVER_H_ */
