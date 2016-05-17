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
#include "Medium.h"
#include "Primitives.h"

namespace tethex {
class Mesh;
} /* namespace tethex */

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;

/** Class to calculate electric field from electric potential */
class LaplacePostProcessor : public DataPostprocessorVector<DIM> {
public:
    LaplacePostProcessor(const string data_name) :
        DataPostprocessorVector<DIM>(data_name, update_values | update_gradients) {
    }

    void compute_derived_quantities_scalar (
            const vector<double>           &uh,
            const vector<Tensor<1,DIM>>    &duh,
            const vector<Tensor<2,DIM>>    &dduh,
            const vector<Point<DIM>>       &normals,
            const vector<Point<DIM>>       &evaluation_points,
            vector<Vector<double>>         &computed_quantities) const {

        Assert(computed_quantities.size() == uh.size(),
                ExcDimensionMismatch(computed_quantities.size(), uh.size()));
        for (unsigned int i = 0; i < computed_quantities.size(); ++i)
            for (unsigned int coord = 0; coord < 2; ++coord)
                computed_quantities[i](coord) = duh[i][coord];
    }
};

/**
 * Class to solve differential equations with finite element method taking as input the FEM mesh
 */
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

    void extract_solution_at_medium(Medium &surf);
    void extract_elfield_at_surf_old(Medium &surf, const string file_name);

    struct Solution {
        vector<Point3d> point;
        vector<Vec3d> elfield;
        vector<double> elfield_norm;
        vector<double> potential;
    };

    Solution solution;

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
    Vector<double> laplace_solution;
    Vector<double> system_rhs;
    ConstraintMatrix constraints;

    const string get_file_type(const string file_name);
    bool on_boundary(const double face, const double face_min, const double face_max);

    double get_potential_at_node(const int &cell_indx, const int &vert_indx);
    double get_potential_at_point(Point<DIM> &point);

    Tensor<1,DIM> get_elfield_at_node(const int &cell, const int &vert_indx);
    Tensor<1,DIM> get_elfield_at_point(Point<DIM> &point);

    vector<int> get_medium2node_map(Medium &medium);
    vector<int> get_node2elem_map();
    vector<int> get_node2vert_map();

    void reserve_solution(const int n_nodes);
    void write_xyz(const string file_name);
};

} /* namespace femocs */

#endif /* DEALII_H_ */
