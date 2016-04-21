/*
 * DealII.h
 *
 *  Created on: 11.2.2016
 *      Author: veske
 */

#ifndef DEALII_H_
#define DEALII_H_

#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include "Macros.h"
#include "AtomReader.h"
#include "Tethex.h"

namespace tethex {
class Mesh;
} /* namespace tethex */

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;

class DealII {
public:
    DealII(const int poly_degree, const double neumann);
    void run();
    void import_file(const string file_name);

    void make_simple_mesh();
    void import_tetgen_mesh(femocs::Mesh* mesh);
    void import_tethex_mesh(tethex::Mesh* mesh);

    void output_mesh(const string file_name);
    void output_results(const string file_name);

    void setup_system();
    void mark_boundary(const AtomReader::Sizes* sizes, const AtomReader::Types* types);
    void assemble_system(const AtomReader::Types* types);
    void solve_umfpack();
    void solve_cg();

private:

    const string get_file_type(const string file_name);
    bool on_boundary(const double face, const double face_min, const double face_max);

    Triangulation<3> triangulation;
    FE_Q<DIM> fe;
    DoFHandler<3> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> solution;
    Vector<double> system_rhs;
    ConstraintMatrix constraints;

    double neumann; //!< gradient length for the solution at some boundary (in current case at the top)
};

} /* namespace femocs */

#endif /* DEALII_H_ */
