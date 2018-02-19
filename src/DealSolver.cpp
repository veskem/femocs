/*
 * DealSolver.cpp
 *
 *  Created on: 12.2.2018
 *      Author: veske
 */

#include "DealSolver.h"

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>


using namespace dealii;
using namespace std;
namespace fch {

template<int dim>
DealSolver<dim>::DealSolver() : fe(shape_degree), dof_handler(triangulation) {}

template<int dim>
DealSolver<dim>::DealSolver(Triangulation<dim> *tria) : fe(shape_degree), dof_handler(*tria) {}

template<int dim>
unsigned int DealSolver<dim>::solve_cg(int max_iter, double tol, double ssor_param) {
    SolverControl solver_control(max_iter, tol);
    SolverCG<> solver(solver_control);

    if (ssor_param > 0.0) {
        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(system_matrix, ssor_param);
        solver.solve(system_matrix, solution, system_rhs, preconditioner);
    } else
        solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    
    solution_old = solution;
    return solver_control.last_step();
}

template<int dim>
vector<double> DealSolver<dim>::shape_funs(const Point<dim> &p, int cell_index) const {
    return shape_funs(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
vector<double> DealSolver<dim>::shape_funs(const Point<dim> &p, const int cell_index,
                                               Mapping<dim,dim>& mapping) const {

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, max(0,cell_index), &dof_handler);
    //point in transformed unit cell coordinates
    Point<dim> p_cell;

    if (cell_index < 0){ // in case the cell index is unknown (argument cell_index < 0)
        const std::pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    } else // cell index is known
        p_cell = mapping.transform_real_to_unit_cell(cell, p);

    Point<dim> p_unit_cell = GeometryInfo<dim>::project_to_unit_cell(p_cell);

    //create virtual quadrature point
    const Quadrature<dim> quadrature(p_unit_cell);

    //define fevalues object
    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);

    fe_values.reinit(cell);

    vector<double> sfuns(fe.dofs_per_cell);

    for (int i = 0; i < sfuns.size(); i++)
        sfuns[i] = fe_values.shape_value(i,0);

    return sfuns;
}

template<int dim>
double DealSolver<dim>::probe_solution(const Point<dim> &p) const {
    return VectorTools::point_value(dof_handler, solution, p);
}

template<int dim>
double DealSolver<dim>::max_solution() const {
    return solution.linfty_norm();
}

template<int dim>
double DealSolver<dim>::get_cell_vol(const int i) const {
    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, i, &dof_handler);
    return cell->measure();
}

template<int dim>
void DealSolver<dim>::import_mesh_from_file(const string &file_name) {
    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.import_mesh_from_file(&triangulation, file_name);
    mark_boundary();
}

template<int dim>
bool DealSolver<dim>::import_mesh_directly(vector<Point<dim>> vertices, vector<CellData<dim>> cells) {
    try {
        SubCellData subcelldata;
        // Do some clean-up on vertices...
        GridTools::delete_unused_vertices(vertices, cells, subcelldata);
        // ... and on cells
        GridReordering<dim, dim>::invert_all_cells_of_negative_grid(vertices, cells);
        // Clean previous mesh
        triangulation.clear();
        // Create new mesh
        triangulation.create_triangulation_compatibility(vertices, cells, SubCellData());
    } catch (exception &exc) {
        return false;
    }

    mark_boundary();
    return true;
}

template<int dim>
void DealSolver<dim>::output_mesh(const string &file_name) {
    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.output_mesh(&triangulation, file_name);
}

template<int dim>
void DealSolver<dim>::get_surface_nodes(vector<Point<dim>>& nodes) const {
    const int n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;
    nodes.clear();

    // Loop over copper interface cells
    typename DoFHandler<dim>::active_cell_iterator cell;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (int f = 0; f < n_faces_per_cell; ++f)
            if (cell->face(f)->boundary_id() == BoundaryId::copper_surface)
                nodes.push_back(cell->face(f)->center());
}

template<int dim>
void DealSolver<dim>::assemble_rhs(const BoundaryId bid) {

    QGauss<dim-1> face_quadrature_formula(this->quadrature_degree);
    FEFaceValues<dim> fe_face_values(this->fe, face_quadrature_formula,
            update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell;

    // Iterate over all cells (quadrangles in 2D, hexahedra in 3D) of the mesh
    unsigned int boundary_face_index = 0;
    for (cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell) {
        cell_rhs = 0;

        // Apply boundary condition at faces on top of vacuum domain
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == bid) {
                fe_face_values.reinit(cell, f);
                double bc_value = get_face_bc(boundary_face_index++);

                // Compose local rhs update
                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += fe_face_values.shape_value(i, q)
                                * bc_value * fe_face_values.JxW(q);
                    }
                }
            }
        }

        // Add the current cell rhs entries to the system rhs
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            this->system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}

template<int dim>
void DealSolver<dim>::assemble_finalize(const BoundaryId bid, const double boundary_val) {
    map<types::global_dof_index, double> boundary_values;
    if (boundary_val > 0.0)
        VectorTools::interpolate_boundary_values(this->dof_handler, bid, ConstantFunction<dim>(boundary_val), boundary_values);
    else
        VectorTools::interpolate_boundary_values(this->dof_handler, bid, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, this->system_matrix, this->solution, this->system_rhs);
}

template class DealSolver<2>;
template class DealSolver<3>;

} /* namespace femocs */
