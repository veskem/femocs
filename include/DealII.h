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

const int DIM = 3;          //!< dimensionality of solver
const int POLY_DEGREE = 1;  //!< polynomial degree of the finite elements (1-linear, 2-quadratic, ...)

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
    DealII();

    const void set_neumann(const double neumann);
    const int get_n_dofs();

    const void run();
    const void import_file(const string file_name);

    const void make_simple_mesh();
    const void import_tetgen_mesh(femocs::Mesh* mesh);
    const void import_tethex_mesh(tethex::Mesh* mesh);

    const void output_mesh(const string file_name);
    const void output_results(const string file_name);

    const void setup_system();
    const void mark_boundary(const AtomReader::Sizes* sizes, const AtomReader::Types* types);
    const void assemble_system(const AtomReader::Types* types);
    const void solve_umfpack();
    const void solve_cg();

    const void extract_solution_at_medium(Medium &surf);
    const void extract_elfield_at_surf_old(Medium &surf, const string file_name);

    const void export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);

    struct Solution {
        vector<int> id;
        vector<Point3> point;
        vector<Vec3> elfield;
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

    const double get_potential_at_node(const int &cell_indx, const int &vert_indx);
    const double get_potential_at_point(Point<DIM> &point);

    const Tensor<1,DIM> get_elfield_at_node(const int &cell, const int &vert_indx);
    const Tensor<1,DIM> get_elfield_at_point(Point<DIM> &point);

    const vector<int> get_medium2node_map(Medium &medium);
    const vector<int> get_node2elem_map();
    const vector<int> get_node2vert_map();

    const void reserve_solution(const int n_nodes);
    const void write_xyz(const string file_name);
};

} /* namespace femocs */

#endif /* DEALII_H_ */
