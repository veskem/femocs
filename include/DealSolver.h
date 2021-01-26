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
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include "Globals.h"
#include "Medium.h"
#include "FileWriter.h"

using namespace dealii;
using namespace std;

namespace femocs {

template<int dim> class PoissonSolver;
template<int dim> class CurrentHeatSolver;
template<int dim> class CurrentSolver;
template<int dim> class HeatSolver;

/** @brief General class to implement FEM solver in Deal.II
 */
template<int dim>
class DealSolver: public FileWriter {
public:

    DealSolver();
    DealSolver(Triangulation<dim> *tria);
    virtual ~DealSolver() {}
    
    /** Provide triangulation object to get access to the mesh data */
    Triangulation<dim>* get_triangulation() { return &triangulation; }

    /** Provide dof_handler object to get access to the mesh data */
    DoFHandler<dim>* get_dof_handler() { return &dof_handler; }

    /** Provide access to solution vector */
    Vector<double>* get_solution() { return &solution; };

    /** Obtain solution values on mesh vertices */
    void export_solution(vector<double> &solution) const;

    /** Obtain MINUS gradient of solution values on mesh vertices */
    void export_solution_grad(vector<Tensor<1, dim>> &grads) const;

    /** Get the solution value at the specified point. NB: Slow! */
    double probe_solution(const Point<dim> &p) const;

    /** Calculate max value from solution vector */
    double max_solution() const;

    /** Check whether the solution is within limits */
    bool check_limits(const double low_limit, const double high_limit);

    /** Values of the shape functions for a point with respect to the cell.
     *  The cell that surrounds the point has to be found before hand. */
    vector<double> shape_funs(const Point<dim> &point, int cell_index) const;

    /** Values of the shape function gradients for a point with respect to the cell.
     *  The cell that surrounds the point has to be found before hand. */
    vector<Tensor<1, dim, double>> shape_fun_grads(const Point<dim> &point, const int cell_index) const;

    /** Return the volume/area of i-th cell */
    double get_cell_vol(const int i) const;

    /** Return number of active cells */
    int get_n_cells() const { return triangulation.n_active_cells(); }

    /** Export mesh vertices into Medium */
    void export_vertices(Medium& medium);

    /** Export mesh dofs into vector */
    void export_dofs(vector<Point<dim>>& points) const;

    /** Export the centroids of surface faces in the order required by assemble_rhs */
    void export_surface_centroids(Medium& medium) const;

    /** Modify DOF solution with nodal solution */
    void import_solution(const vector<double>* new_solution);

    /** Return # degrees of freedom of system */
    int size() const { return dof_handler.n_dofs(); }

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

    /** Print the statistics about the mesh and # degrees of freedom */
    friend ostream& operator <<(ostream &os, const DealSolver<dim>& d) {
        os << "#elems=" << d.tria->n_active_cells()
                << ", #faces=" << d.tria->n_active_faces()
                << ", #edges=" << d.tria->n_active_lines()
                << ", #nodes=" << d.tria->n_used_vertices()
                << ", #dofs=" << d.dof_handler.n_dofs();
        return os;
    }

    /** Return mesh statistics as string */
    string to_str() const { ostringstream ss; ss << (*this); return ss.str(); }

    struct Stat {
        double sol_min;   ///< minimum solution value
        double sol_max;   ///< maximum solution value

        friend ostream& operator <<(ostream &os, const DealSolver<dim>::Stat &s) {
            os << "min=" << s.sol_min << ", max=" << s.sol_max;
            return os;
        }
    } stat;

protected:
    static constexpr unsigned int shape_degree = 1;   ///< degree of the shape functions (linear, quadratic etc elements)
    static constexpr unsigned int quadrature_degree = shape_degree + 1;  ///< degree of the Gaussian numerical integration

    double dirichlet_bc_value;

    Triangulation<dim> triangulation;        ///< local triangulation
    Triangulation<dim>* tria;                ///< pointer to external triangulation
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern sparsity_pattern;        ///< structure for sparse matrix representation
    SparseMatrix<double> system_matrix;      ///< system matrix of matrix equation
    Vector<double> system_rhs;               ///< right-hand-side of the matrix equation
    Vector<double> solution;                 ///< resulting solution in the mesh nodes

    vector<double> dof_volume;               ///< integral of the shape functions
    vector<unsigned> vertex2dof;             ///< map of vertex to dof indices
    vector<unsigned> vertex2cell;            ///< map of vertex to cell indices
    vector<unsigned> vertex2node;            ///< map of vertex to cell node indices

    /** Variables used during the assembly of matrix equation */
    map<types::global_dof_index, double> boundary_values;

    /** Data holding copy of global matrix & rhs for parallel assembly */
    struct LinearSystem {
        Vector<double>* global_rhs;
        SparseMatrix<double>* global_matrix;
        LinearSystem(Vector<double>* rhs, SparseMatrix<double>* matrix);
    };

    /** Data for parallel local matrix & rhs assembly */
    struct ScratchData {
        FEValues<dim> fe_values;
        ScratchData(const FiniteElement<dim> &fe, const Quadrature<dim> &quadrature, const UpdateFlags flags);
        ScratchData(const ScratchData &scratch_data);
    };

    /** Data for coping local matrix & rhs into global one during parallel assembly */
    struct CopyData {
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;
        vector<unsigned int> dof_indices;
        unsigned int n_dofs, n_q_points;
        CopyData(const unsigned dofs_per_cell, const unsigned n_q_points);
    };

    /** Copy the matrix & rhs vector contribution of a cell into global matrix & rhs vector */
    // Only one instance of this function should be running at a time!
    void copy_global_cell(const CopyData &copy_data, LinearSystem &system) const;

    /** Helper function for the public shape_funs */
    vector<double> shape_funs(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Helper function for the public shape_fun_grads */
    vector<Tensor<1, dim, double>> shape_fun_grads(const Point<dim> &point, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Solve the matrix equation using Conjugate Gradient method
     * @param n_steps     maximum number of iterations allowed
     * @param tolerance   tolerance of the solution
     * @param ssor_param  parameter to SSOR preconditioner. Its fine tuning optimises calculation time
     */
    int solve_cg(const int n_steps, const double tolerance, const double ssor_param);

    /** Set up dynamic sparsity pattern for calculations */
    void setup_system();

    /** Modify the right-hand-side vector of the matrix equation */
    void assemble_rhs(const int bid);

    /** Give the value to all DOFs with given boundary ID */
    void append_dirichlet(const int bid, const double value);

    /** Apply all dirichlet boundary conditions to the system.
     * This should be the last function call to setup the equations, before calling solve(). */
    void apply_dirichlet();

    /** Calculate the volumes (dim=3) or areas (dim=2) of dofs used during integration */
    void calc_dof_volumes();

    /** Calculate mapping between vertex and dof indices */
    void calc_vertex2dof();

    /** Specify file types that can be written */
    bool valid_extension(const string &ext) const {
        return ext == "vtk" || ext == "vtks" || ext == "msh";
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
    void mark_boundary(int top, int bottom, int sides, int other);

    friend class PoissonSolver<dim>;
    friend class CurrentHeatSolver<dim>;
    friend class CurrentSolver<dim>;
    friend class HeatSolver<dim>;
};

} /* namespace femocs */

#endif /* DEALSOLVER_H_ */
