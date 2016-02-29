/*
 * DealII.h
 *
 *  Created on: 11.2.2016
 *      Author: veske
 */

#ifndef DEALII_H_
#define DEALII_H_

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/base/timer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/grid/grid_reordering.h>

#include <memory>
#include "Femocs.h"
#include "Mesh.h"
#include "tethex.h"

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;

class DealII {
public:
    DealII(const int poly_degree, const double neumann, Femocs::SimuCell* simucell);
    void run();
    void import_file(const string file_name);
    void make_simple_mesh();
    void import_tetgen_mesh(shared_ptr<Mesh> mesh);
    void import_tethex_mesh(tethex::Mesh* mesh);
    void import_tethex_mesh_old(tethex::Mesh* mesh);

    void output_mesh(const string file_name);
    void output_results(const string file_name);

    void setup_system();
    void mark_boundary();
    void assemble_system();
    void solve_umfpack();
    void solve_cg();
    void distort_solution(const double dist_ampl);
    void distort_solution_const(const double dist_ampl);
    void distort_solution_one(const double dist_ampl, const int i);


private:

    const string get_file_type(const string& file_name);
    bool on_boundary(const double face, const double face_min, const double face_max);

    Triangulation<3> triangulation;
    FE_Q<DIM> fe;
    DoFHandler<3> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> solution;
    Vector<double> system_rhs;
    ConstraintMatrix constraints;

    //Declare enumerators for domain distinction. To be implemented when we apply multiphysics.
    struct Types {
        unsigned int vacuum;
        unsigned int bulk;
        unsigned int top;
        unsigned int surface;
        unsigned int side;
    };

    struct Sizes {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        double zmaxbox;
        double zminbox;
    };

    Types types;
    Sizes sizes;
    double neumann;  //!< gradient length for the solution at some boundary (in current case at the top)
};

} /* namespace femocs */

#endif /* DEALII_H_ */
