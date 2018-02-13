/*
 * DealSolver.cpp
 *
 *  Created on: 12.2.2018
 *      Author: veske
 */

#include "DealSolver.h"


using namespace dealii;
using namespace std;
namespace fch {

template<int dim>
DealSolver<dim>::DealSolver() : fe(shape_degree), dof_handler(triangulation) {}

template<int dim>
DealSolver<dim>::DealSolver(Triangulation<dim> *tria_) : fe(shape_degree), dof_handler(*tria_) {}

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
void DealSolver<dim>::get_surface_nodes(vector<Point<dim>>& nodes) {
    const int n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;
    nodes.clear();

    // Loop over copper interface cells
    typename DoFHandler<dim>::active_cell_iterator cell;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (int f = 0; f < n_faces_per_cell; ++f)
            if (cell->face(f)->boundary_id() == BoundaryId::copper_surface)
                nodes.push_back(cell->face(f)->center());
}

template class DealSolver<2>;
template class DealSolver<3>;

} /* namespace femocs */
