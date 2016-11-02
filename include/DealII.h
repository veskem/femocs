/*
 * DealII.h
 *
 *  Created on: 11.2.2016
 *      Author: veske, Kristjan Eimre, Robert Aare
 */

#ifndef DEALII_H_
#define DEALII_H_

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

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;         //!< dimensionality of solver
const int POLY_DEGREE = 1; //!< polynomial degree of the finite elements (1-linear, 2-quadratic etc)

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
            for (unsigned int coord = 0; coord < DIM; ++coord)
                computed_quantities[i](coord) = duh[i][coord];
    }
};

/** Class to solve differential equations with finite element method taking as input the FEM mesh */
class DealII {
public:
    DealII();

    /** Specify the Neumann boundary condition value */
    const void set_neumann(const double neumann);

    /** Import mesh from file */
    const bool import_mesh(const string &file_name);

    /** Import vertices, quadrangles and hexahedra */
    const bool import_mesh(tethex::Mesh& mesh);

    /** Import vertices and hexahedra and ignore quadrangles */
    const bool import_mesh_wo_faces(tethex::Mesh& mesh);

    /** Make mesh 4x denser in the sphere with the centre in origin and radius r_cut */
    const void refine_mesh(const Point3 &origin, const double r_cut);

    /** Write the mesh to file */
    const void write_mesh(const string &file_name);

    /** Write the calculated electric potential and electric field to file */
    const void write_results(const string &file_name);

    /** Mark boundary faces, distribute degrees of freedom  and initialise data */
    const void setup_system(const AtomReader::Sizes& sizes);

    /** Insert boundary conditions to the system */
    const void assemble_system();

    /** Run the calculation with UMFPACK solver */
    const void solve_umfpack();

    /** Run the calculation with Conjugate Gradient solver */
    const void solve_cg();

    /** Calculate electric potential at an arbitrary point inside the mesh.
     * Point outside the mesh gives an error. */
    const double get_potential(const double x, const double y, const double z);

    /** Calculate electric potential at a set of mesh nodes
     * * @param cell_indxs - indices of elements where the node is located
     * @param vert_indxs - indices of vertices of the element where the node is located */
    const vector<double> get_potential(const vector<int> &cell_indxs, const vector<int> &vert_indxs);

    /** Calculate electric field at an arbitrary point inside the mesh.
     * Point outside the mesh gives an error. */
    const Vec3 get_elfield(const double x, const double y, const double z);

    /** Calculate electric field at a set of mesh nodes
     * @param cell_indxs - indices of elements where the node is located
     * @param vert_indxs - indices of vertices of the element where the node is located */
    const vector<Vec3> get_elfield(const vector<int> &cell_indxs, const vector<int> &vert_indxs);

    /** string stream prints the statistics about the system */
    friend std::ostream& operator <<(std::ostream &os, const DealII& d) {
        os << "#elems=" << d.triangulation.n_active_cells()
                << ",\t#faces=" << d.triangulation.n_active_faces()
                << ",\t#edges=" << d.triangulation.n_active_lines()
                << ",\t#nodes=" << d.triangulation.n_used_vertices()
                << ",\t#dofs=" << d.dof_handler.n_dofs();
        return os;
    }

    const unsigned int n_verts_per_elem = GeometryInfo<3>::vertices_per_cell; //!< # vertices in element
    const unsigned int n_verts_per_face = GeometryInfo<2>::vertices_per_cell; //!< # vertices in face

    Triangulation<DIM> triangulation;
    DoFHandler<DIM> dof_handler;

private:
    double neumann;         //!< Neumann boundary condition value applied to upper region faces

    FE_Q<DIM> fe;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> laplace_solution;
    Vector<double> system_rhs;
    
    /** Mark the boundary faces of the mesh */
    const void mark_boundary_faces(const AtomReader::Sizes& sizes);
};

} /* namespace femocs */

#endif /* DEALII_H_ */
