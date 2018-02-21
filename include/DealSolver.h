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

#include <fstream>
#include "Globals.h"

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
    
    /** Provide triangulation object to get access to the mesh data */
    Triangulation<dim>* get_triangulation() { return &triangulation; }

    /** Provide dof_handler object to get access to the mesh data */
    DoFHandler<dim>* get_dof_handler() { return &dof_handler; }

    /** Get the solution value at the specified point. NB: Slow! */
    double probe_solution(const Point<dim> &p) const;

    /** Calculate max value from solution vector */
    double max_solution() const;

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
    bool import_mesh(const string &file_name);

    /**
     * imports mesh directly from vertex and cell data and sets the boundary indicators
     * @return true if success, otherwise false
     */
    bool import_mesh(vector<Point<dim>> vertices, vector<CellData<dim>> cells);

    /** Write data to file; data type and format is determined with the file extention */
    void write(const string &filename) const;

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

    double dirichlet_bc_value;

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;        ///< structure for sparse matrix representation
    SparseMatrix<double> system_matrix;      ///< system matrix of matrix equation
    SparseMatrix<double> system_matrix_save; ///< system matrix of matrix equation (save before Dirichlet BCs remove dofs)
    Vector<double> system_rhs;            ///< right-hand-side of the matrix equation

    Vector<double> solution;              ///< resulting solution in the mesh nodes
    Vector<double> solution_save;         ///< resulting solution in the mesh nodes

    /** Variables used during the assembly of matrix equation */
    map<types::global_dof_index, double> boundary_values;

    /** Helper function for the public shape_funs */
    vector<double> shape_funs(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Solve the matrix equation using Conjugate Gradient method
     * @param n_steps     maximum number of iterations allowed
     * @param tolerance   tolerance of the solution
     * @param ssor_param  parameter to SSOR preconditioner. Its fine tuning optimises calculation time
     */
    unsigned int solve_cg(const int n_steps, const double tolerance, const double ssor_param);

    /** Set up dynamic sparsity pattern for calculations */
    void setup_system();

    /** Modify the right-hand-side vector of the matrix equation */
    void assemble_rhs(const BoundaryId bid);

    /** Give the value to all DOFs with given boundary ID */
    void append_dirichlet(const BoundaryId bid, const double value);

    /** Apply all dirichlet boundary conditions to the system.
     * This should be the last function call to setup the equations, before calling solve(). */
    void apply_dirichlet();

    /** saves the system matrix (useful before BCs have been applied) */
    void save_system() {
        system_matrix_save.copy_from(system_matrix);
    }

    /** Restores the system matrix to the saved one */
    void restore_system() {
        system_matrix.copy_from(system_matrix_save);
    }

    /** Write the mesh to msh file that can in turn be imported to Deal.II */
    void write_msh(ofstream& out) const;

    /** Write the simulation results to file in vtk format that can be visualized in Paraview */
    virtual void write_vtk(ofstream& out) const;

    /** Mark different regions of the mesh */
    virtual void mark_mesh() {}

    /** Return the boundary condition value at the centroid of face */
    virtual double get_face_bc(const unsigned int face) const { return 1.0; }

    /**
     * Marks boundary (top, bottom, sides, vacuum-copper surface) faces (3D) or edges (2D) in the mesh
     * @param top     id for top faces|edges
     * @param bottom  id for bottom faces|edges
     * @param sides   id for lateral faces|edges
     * @param other   id for vacuum-copper surface faces|edges
     */
    void mark_boundary(char top, char bottom, char sides, char other);

    friend class PoissonSolver<dim>;
    friend class CurrentHeatSolver<dim>;
    friend class CurrentSolver<dim>;
    friend class HeatSolver<dim>;
};

} /* namespace femocs */

#endif /* DEALSOLVER_H_ */
