/*
 * DealII.h
 *
 *  Created on: 11.2.2016
 *      Author: veske, Kristjan Eimre, Robert Aare
 */

#ifndef DEALII_H_
#define DEALII_H_

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/data_out.h>

#include "Macros.h"
#include "Medium.h"

using namespace std;
using namespace dealii;
namespace femocs {

const int DIM = 3;         ///< dimensionality of solver
const int POLY_DEGREE = 1; ///< polynomial degree of the finite elements (1-linear, 2-quadratic etc)

/** Class to calculate electric field from electric potential */
class LaplacePostProcessor : public DataPostprocessorVector<DIM> {
public:
    LaplacePostProcessor(const string &data_name) :
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
    void set_applied_efield(const double neumann);

    /** Import mesh from file */
    bool import_mesh(const string &file_name);

    /** Import vertices and hexahedra into Deal.II */
    bool import_mesh(vector<Point<DIM>> vertices, vector<CellData<DIM>> cells);

    /** Make mesh 4x denser in the sphere */
    void refine_mesh(const Point<DIM> &origin, const double radius);

    /** Write the calculated electric potential and electric field to file */
    void write(const string &file_name);

    /** Write the mesh to file */
    void write_mesh(const string &file_name);

    /** Mark boundary faces, distribute degrees of freedom  and initialise data */
    void setup_system(const Medium::Sizes& sizes);


    Triangulation<DIM>* get_triangulation();

    DoFHandler<DIM>* get_dof_handler();


    /** Insert boundary conditions to the system */
    void assemble_system();

    /** Run the calculation with UMFPACK solver */
    void solve_umfpack();

    /** Run the calculation with Conjugate Gradient solver */
    void solve_cg();

    /** Calculate electric potential at an arbitrary point inside the mesh.
     * Point outside the mesh gives an error. */
    double get_potential(const double x, const double y, const double z);

    /** Calculate electric potential at a set of mesh nodes
     * * @param cell_indxs - indices of elements where the node is located
     * @param vert_indxs - indices of vertices of the element where the node is located */
    vector<double> get_potential(const vector<int> &cell_indxs, const vector<int> &vert_indxs);

    /** Calculate electric field at an arbitrary point inside the mesh.
     * Point outside the mesh gives an error. */
    Vec3 get_efield(const double x, const double y, const double z);

    /** Calculate electric field at a set of mesh nodes
     * @param cell_indxs - indices of elements where the node is located
     * @param vert_indxs - indices of vertices of the element where the node is located */
    vector<Tensor<1, DIM>> get_efield(const vector<int> &cell_indxs, const vector<int> &vert_indxs);

    /** string stream prints the statistics about the system */
    friend std::ostream& operator <<(std::ostream &os, const DealII& d) {
        os << "#elems=" << d.triangulation.n_active_cells()
                << ",\t#faces=" << d.triangulation.n_active_faces()
                << ",\t#edges=" << d.triangulation.n_active_lines()
                << ",\t#nodes=" << d.triangulation.n_used_vertices()
                << ",\t#dofs=" << d.dof_handler.n_dofs();
        return os;
    }

    const unsigned int n_verts_per_elem = GeometryInfo<3>::vertices_per_cell; ///< # vertices in element
    const unsigned int n_verts_per_face = GeometryInfo<2>::vertices_per_cell; ///< # vertices in face

    Triangulation<DIM> triangulation;
    DoFHandler<DIM> dof_handler;

private:
    double neumann;         ///< Neumann boundary condition value applied to upper region faces

    FE_Q<DIM> fe;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> laplace_solution;
    Vector<double> system_rhs;
    
    /** Mark the boundary faces of the mesh */
    void mark_boundary_faces(const Medium::Sizes& sizes);
};

} /* namespace femocs */

#endif /* DEALII_H_ */
