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
#include "Media.h"

namespace tethex {
class Mesh;
} /* namespace tethex */

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;

class LaplacePostProcessor : public DataPostprocessorVector<DIM> {
public:
    LaplacePostProcessor(const string data_name);

    virtual void compute_derived_quantities_scalar (
            const vector<double>           &uh,
            const vector<Tensor<1,DIM>>    &duh,
            const vector<Tensor<2,DIM>>    &dduh,
            const vector<Point<DIM>>       &normals,
            const vector<Point<DIM>>       &evaluation_points,
            vector<Vector<double>>         &computed_quantities) const;
};

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

    void probe_results(const int N);
    void extract_elfield_at_surf(Surface* surf, const string file_name);

private:
    const unsigned int n_verts_per_elem = GeometryInfo<DIM>::vertices_per_cell;
    const unsigned int n_verts_per_face = GeometryInfo<DIM-1>::vertices_per_cell;
    const unsigned int n_faces_per_elem = GeometryInfo<DIM>::faces_per_cell;

    double neumann;         //!< Gradient length for the solution at top boundary

    Triangulation<3> triangulation;
    FE_Q<DIM> fe;
    DoFHandler<3> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> potential;
    Vector<double> system_rhs;
    ConstraintMatrix constraints;

    std::vector<Tensor<1,DIM>> elfield;
    std::vector<double> elfield_norm;

    const string get_file_type(const string file_name);
    bool on_boundary(const double face, const double face_min, const double face_max);

    double get_potential_at_node(Point<DIM> &node);
    double get_potential_at_point(Point<DIM> &point);

    Tensor<1,DIM> get_elfield_at_node(Point<DIM> &node, const unsigned int &cell, const unsigned int &vert_indx);
    Tensor<1,DIM> get_elfield_at_node_old(Point<DIM> &node);

    Tensor<1,DIM> get_elfield_at_point(Point<DIM> &point);

    vector<unsigned int> get_node2elem_map(const int n_surf);
};

} /* namespace femocs */

#endif /* DEALII_H_ */
