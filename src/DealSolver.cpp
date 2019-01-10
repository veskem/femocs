/*
 * DealSolver.cpp
 *
 *  Created on: 12.2.2018
 *      Author: veske
 */

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_tools.h>

#include "DealSolver.h"
#include "Macros.h"
#include "Globals.h"


using namespace dealii;
using namespace std;

namespace femocs {

template<int dim>
DealSolver<dim>::DealSolver() :
        dirichlet_bc_value(0), tria(&triangulation), fe(shape_degree), dof_handler(triangulation) {}

template<int dim>
DealSolver<dim>::DealSolver(Triangulation<dim> *tr) :
        dirichlet_bc_value(0), tria(tr), fe(shape_degree), dof_handler(*tr) {}

template<int dim>
DealSolver<dim>::LinearSystem::LinearSystem(Vector<double>* rhs, SparseMatrix<double>* matrix) :
    global_rhs(rhs), global_matrix(matrix)
{}

template<int dim>
DealSolver<dim>::ScratchData::ScratchData (const FiniteElement<dim>  &fe,
        const Quadrature<dim> &quadrature, const UpdateFlags ul) :
    fe_values(fe, quadrature, ul)
{}

template<int dim>
DealSolver<dim>::ScratchData::ScratchData (const ScratchData &sd):
    fe_values(sd.fe_values.get_fe(), sd.fe_values.get_quadrature(), sd.fe_values.get_update_flags())
{}

template<int dim>
DealSolver<dim>::CopyData::CopyData(const unsigned dofs_per_cell, const unsigned n_qp):
    cell_matrix(dofs_per_cell, dofs_per_cell),
    cell_rhs(dofs_per_cell),
    dof_indices(dofs_per_cell),
    n_dofs(dofs_per_cell), n_q_points(n_qp)
{}

template<int dim>
void DealSolver<dim>::copy_global_cell(const CopyData &copy_data, LinearSystem &system) const {
    system.global_rhs->add(copy_data.dof_indices, copy_data.cell_rhs);

    for (unsigned int i = 0; i < copy_data.n_dofs; ++i) {
        for (unsigned int j = 0; j < copy_data.n_dofs; ++j)
            system.global_matrix->add(copy_data.dof_indices[i], copy_data.dof_indices[j], copy_data.cell_matrix(i, j));
    }
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

    if (cell_index < 0) { // in case the cell index is unknown (argument cell_index < 0)
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

    // store shape functions
    vector<double> sfuns(fe.dofs_per_cell);
    for (unsigned int i = 0; i < sfuns.size(); i++)
        sfuns[i] = fe_values.shape_value(i,0);

    return sfuns;
}

template<int dim>
vector<Tensor<1, dim, double>> DealSolver<dim>::shape_fun_grads(const Point<dim> &p, const int cell_index) const {
    return shape_fun_grads(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
vector<Tensor<1, dim, double>> DealSolver<dim>::shape_fun_grads(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const {

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&this->triangulation, 0, max(0,cell_index), &this->dof_handler);

    // transform the point from real to unit cell coordinates
    Point<dim> p_cell;
    if (cell_index < 0) {
        const pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, this->dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    } else
        p_cell = mapping.transform_real_to_unit_cell(cell, p);

    const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));

    FEValues<dim> fe_values(mapping, this->fe, quadrature, update_gradients);
    fe_values.reinit(cell);

    // store shape function gradients
    vector<Tensor<1, dim, double>> sfun_grads(fe.dofs_per_cell);
    for (unsigned int i = 0; i < sfun_grads.size(); i++)
        sfun_grads[i] = fe_values.shape_grad(i,0);

    return sfun_grads;
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
bool DealSolver<dim>::import_mesh(const string &file_name) {
    const string file_type = get_file_type(file_name);
    require(file_type == "msh", "Unimplemented file type for mesh importing: " + file_type);

    ifstream infile(file_name);
    require(infile.is_open(), "Can't open a file " + file_name);

    GridIn<dim, dim> gi;
    gi.attach_triangulation(triangulation);
    gi.read_msh(infile);

    mark_mesh();
    return true;
}

template<int dim>
bool DealSolver<dim>::import_mesh(vector<Point<dim>> vertices, vector<CellData<dim>> cells) {
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

    mark_mesh();
    return true;
}

template<int dim>
void DealSolver<dim>::write_vtk(ofstream& out) const {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    data_out.write_vtk(out);
}

template<int dim>
void DealSolver<dim>::write_msh(ofstream& out) const {
    GridOut grid_out;
    grid_out.set_flags(GridOutFlags::Msh(true, true));
    grid_out.write_msh(triangulation, out);
}

template<int dim>
void DealSolver<dim>::export_surface_centroids(femocs::Medium& medium) const {
    const int n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;
    typename DoFHandler<dim>::active_cell_iterator cell;

    unsigned int n_nodes = 0;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (int f = 0; f < n_faces_per_cell; ++f)
            if (cell->face(f)->boundary_id() == BoundaryID::copper_surface)
                n_nodes++;

    medium.reserve(n_nodes);

    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (int f = 0; f < n_faces_per_cell; ++f)
            if (cell->face(f)->boundary_id() == BoundaryID::copper_surface)
                medium.append( femocs::Point3(cell->face(f)->center()) );
}

template<int dim>
void DealSolver<dim>::export_vertices(Medium& medium) {
    vector<Point<dim>> support_points;
    export_dofs(support_points);

    const unsigned int n_verts = tria->n_used_vertices();
    this->calc_vertex2dof();

    require(vertex2dof.size() == n_verts, "Mismatch between vertex2dof size and #vertices: "
            + d2s(vertex2dof.size()) + " vs " + d2s(n_verts));

    medium.reserve(n_verts);
    for (int i = 0; i < n_verts; ++i)
        medium.append( Atom(i, support_points[vertex2dof[i]], 0) );
}

template<int dim>
void DealSolver<dim>::export_dofs(vector<Point<dim>>& points) const {
    points.resize(size());
    DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping,
            dof_handler, points);
}

template<int dim>
void DealSolver<dim>::get_nodal_solution(vector<double>& solution) {
    const unsigned int n_verts = tria->n_used_vertices();
    this->calc_vertex2dof();

    require(vertex2dof.size() == n_verts, "Mismatch between vertex2dof size and #vertices: "
            + d2s(vertex2dof.size()) + " vs " + d2s(n_verts));

    solution.resize(n_verts);
    for (unsigned int i = 0; i < n_verts; ++i)
        solution[i] = this->solution[vertex2dof[i]];
}

template<int dim>
void DealSolver<dim>::set_nodal_solution(const vector<double>* new_solution) {
    const unsigned int n_verts = vertex2dof.size();

    require(new_solution, "Can't use NULL solution vector!");
    require(n_verts == new_solution->size(), "Mismatch between #vertices and solution vector size: "
            + d2s(n_verts) + " vs " + d2s(new_solution->size()));

    // Initialize the solution with non-constant values
    for (size_t i = 0; i < n_verts; i++) {
        this->solution[vertex2dof[i]] = (*new_solution)[i];
    }
}

template<int dim>
void DealSolver<dim>::solution_at(vector<double> &sols,
        const vector<int> &cells, const vector<int> &verts) const
{
    const int n_nodes = cells.size();
    require(n_nodes == verts.size(), "Invalid vectors sizes for cells and vertices: "
            + d2s(n_nodes) + ", " + d2s(verts.size()));

    sols.resize(n_nodes);

    for (unsigned i = 0; i < n_nodes; i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(tria, 0, cells[i], &dof_handler);
        sols[i] = solution[dof_cell->vertex_dof_index(verts[i], 0)];
    }
}

template<int dim>
void DealSolver<dim>::solution_grad_at(vector<Tensor<1, dim>> &grads,
        const vector<int> &cells, const vector<int> &verts) const
{
    const int n_nodes = cells.size();
    require(n_nodes == verts.size(), "Invalid vectors sizes for cells and vertices: "
            + d2s(n_nodes) + ", " + d2s(verts.size()));

    QGauss<dim> quadrature_formula(this->quadrature_degree);
    FEValues<dim> fe_values(this->fe, quadrature_formula, update_gradients);

    vector<Tensor<1, dim>> solution_gradients(quadrature_formula.size());
    grads.resize(n_nodes);

    for (unsigned i = 0; i < n_nodes; i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(tria, 0, cells[i], &dof_handler);

        fe_values.reinit(dof_cell);
        fe_values.get_function_gradients(this->solution, solution_gradients);
        grads[i] = -1.0 * solution_gradients.at(verts[i]);
    }
}

template<int dim>
void DealSolver<dim>::calc_vertex2dof() {
    static constexpr int n_verts_per_elem = GeometryInfo<dim>::vertices_per_cell;
    require(tria, "Pointer to triangulation missing!");
    const unsigned int n_verts = tria->n_used_vertices();
    require(n_verts > 0, "Can't generate map with empty triangulation!");

    // create mapping from mesh vertex to cell index & cell node
    vector<unsigned> vertex2hex(n_verts), vertex2node(n_verts);

    typename DoFHandler<dim>::active_cell_iterator cell;
    for (cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell)
        for (int i = 0; i < n_verts_per_elem; ++i) {
            vertex2hex[cell->vertex_index(i)] = cell->active_cell_index();
            vertex2node[cell->vertex_index(i)] = i;
        }

    // create mapping from vertex index to dof index
    this->vertex2dof.resize(n_verts);
    for (unsigned i = 0; i < n_verts; ++i) {
        typename DoFHandler<dim>::active_cell_iterator cell(tria, 0, vertex2hex[i], &this->dof_handler);
        this->vertex2dof[i] = cell->vertex_dof_index(vertex2node[i], 0);
    }
}

template<int dim>
void DealSolver<dim>::setup_system() {
    require(tria->n_used_vertices() > 0, "Can't setup system with no mesh!");

    this->dof_handler.distribute_dofs(this->fe);
    this->boundary_values.clear();

    const unsigned int n_dofs = size();

    DynamicSparsityPattern dsp(n_dofs);
    DoFTools::make_sparsity_pattern(this->dof_handler, dsp);
    this->sparsity_pattern.copy_from(dsp);

    this->system_matrix.reinit(this->sparsity_pattern);
    this->system_rhs.reinit(n_dofs);
    this->solution.reinit(n_dofs);
    this->solution = this->dirichlet_bc_value;
}

template<int dim>
void DealSolver<dim>::assemble_rhs(const int bid) {

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
        // Loop over all faces (lines in 2D, quadrangles in 3D) of the cell
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            // Apply boundary condition at faces on top of vacuum domain
            if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == bid) {
                fe_face_values.reinit(cell, f);
                double bc_value = get_face_bc(boundary_face_index++);

                // Compose local rhs update
                cell_rhs = 0;
                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += fe_face_values.shape_value(i, q)
                                * bc_value * fe_face_values.JxW(q);
                    }
                }

                // Add the current cell rhs entries to the system rhs
                cell->get_dof_indices(local_dof_indices);
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    this->system_rhs(local_dof_indices[i]) += cell_rhs(i);
            }
        }
    }
}

template<int dim>
void DealSolver<dim>::append_dirichlet(const int bid, const double value) {
    VectorTools::interpolate_boundary_values(this->dof_handler, bid, ConstantFunction<dim>(value), boundary_values);
}

template<int dim>
void DealSolver<dim>::apply_dirichlet() {
    MatrixTools::apply_boundary_values(boundary_values, this->system_matrix, this->solution, this->system_rhs);
}

template<int dim>
int DealSolver<dim>::solve_cg(int max_iter, double tol, double ssor_param) {
    SolverControl solver_control(max_iter, tol);
    SolverCG<> solver(solver_control);
    try {
        if (ssor_param > 0.0) {
            PreconditionSSOR<> preconditioner;
            preconditioner.initialize(system_matrix, ssor_param);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
        } else
            solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

        return solver_control.last_step();
    } catch (exception &exc) {
        return -1 * solver_control.last_step();
    }
}

template<int dim>
void DealSolver<dim>::mark_boundary(int top, int bottom, int sides, int other) {
    static constexpr double eps = 1e-6;
    double xmax = -1e16, ymax = -1e16, zmax = -1e16;
    double xmin = 1e16, ymin = 1e16, zmin = 1e16;

    typename Triangulation<dim>::active_face_iterator face;

    // Loop through the faces and find maximum and minimum values for coordinates
    for (face = triangulation.begin_face(); face != triangulation.end_face(); ++face) {
        if (face->at_boundary()) {
            double x = face->center()[0];
            double y = face->center()[1];
            xmax = max(x, xmax);
            xmin = min(x, xmin);
            ymax = max(y, ymax);
            ymin = min(y, ymin);

            if (dim == 3) {
                double z = face->center()[2];
                zmax = max(z, zmax);
                zmin = min(z, zmin);
            }
        }
    }

    // Loop through the faces and mark them by their position
    for (face = triangulation.begin_face(); face != triangulation.end_face(); ++face) {
        if (face->at_boundary()) {
            if (dim == 2) {
                double x = face->center()[0];
                double y = face->center()[1];

                if (on_boundary(x, xmin, xmax, eps))
                    face->set_all_boundary_ids(sides);
                else if (on_boundary(y, ymax, eps))
                    face->set_all_boundary_ids(top);
                else if (on_boundary(y, ymin, eps))
                    face->set_all_boundary_ids(bottom);
                else
                    face->set_all_boundary_ids(other);
            }
            else if (dim == 3) {
                double x = face->center()[0];
                double y = face->center()[1];
                double z = face->center()[2];

                if (on_boundary(x, xmin, xmax, eps) || on_boundary(y, ymin, ymax, eps))
                    face->set_all_boundary_ids(sides);
                else if (on_boundary(z, zmax, eps))
                    face->set_all_boundary_ids(top);
                else if (on_boundary(z, zmin, eps))
                    face->set_all_boundary_ids(bottom);
                else
                    face->set_all_boundary_ids(other);
            }
        }
    }
}

template<int dim>
void DealSolver<dim>::calc_dof_volumes() {

    QGauss<dim> quadrature_formula(quadrature_degree);
    FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points | update_JxW_values);

    // reset volumes
    dof_volume.reinit(size());
    dof_volume = 0;
    vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell;
    // Iterate over all cells
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell) {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        // Iterate through quadrature points to integrate
        for (unsigned q = 0; q < quadrature_formula.size(); ++q) {
            //iterate through local dofs
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                dof_volume[local_dof_indices[i]] +=  fe_values.JxW(q);
        }
    }
}

template class DealSolver<3>;

} /* namespace femocs */
