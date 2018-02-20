/*
 * currents_and_heating_stationary.cc -> CurrentHeatStatSolver.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan
 */

#include "CurrentsAndHeatingStationary.h"
#include "Utility.h"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>	// UMFpack

#include <cassert>
#include <algorithm>

using namespace dealii;
using namespace std;
namespace fch {

// ----------------------------------------------------------------------------------------
// Class for outputting the current density distribution (calculated from potential distr.)
template<int dim>
class CurrentPostProcessorStat: public DataPostprocessorVector<dim> {
    PhysicalQuantities *pq;
public:
    CurrentPostProcessorStat(PhysicalQuantities *pq_) :
            DataPostprocessorVector<dim>("current_density", update_values | update_gradients), pq(
                    pq_) {
    }

    void compute_derived_quantities_vector(const vector<Vector<double> > & uh,
            const vector<vector<Tensor<1, dim> > > & duh,
            const vector<vector<Tensor<2, dim> > > & /*dduh*/,
            const vector<Point<dim> > & /*normals*/,
            const vector<Point<dim> > & /*evaluation_points*/,
            vector<Vector<double> > & computed_quantities) const {
        for (unsigned int i = 0; i < computed_quantities.size(); i++) {
            double t = uh[i][1]; // temperature
            double sigma = pq->sigma(t);
            for (unsigned int d = 0; d < dim; ++d) {
                double e_field = -duh[i][0][d]; // gradient of the 0-th vector (i.e. potential)
                computed_quantities[i](d) = sigma * e_field;
            }
        }
    }
};

// Class for outputting the electrical conductivity distribution
template<int dim>
class SigmaPostProcessorStat: public DataPostprocessorScalar<dim> {
    PhysicalQuantities *pq;
public:
    SigmaPostProcessorStat(PhysicalQuantities *pq_) :
            DataPostprocessorScalar<dim>("sigma", update_values), pq(pq_) {
    }

    void compute_derived_quantities_vector(const vector<Vector<double> > & uh,
            const vector<vector<Tensor<1, dim> > > & /*duh*/,
            const vector<vector<Tensor<2, dim> > > & /*dduh*/,
            const vector<Point<dim> > & /*normals*/,
            const vector<Point<dim> > & /*evaluation_points*/,
            vector<Vector<double> > & computed_quantities) const {
        for (unsigned int i = 0; i < computed_quantities.size(); i++) {
            double t = uh[i][1]; // temperature
            computed_quantities[i](0) = pq->sigma(t);
        }
    }
};
// ----------------------------------------------------------------------------------------

template<int dim>
CurrentsAndHeatingStationary<dim>::CurrentsAndHeatingStationary() :
        ambient_temperature(ambient_temperature_default), fe(FE_Q<dim>(currents_degree), 1, // Finite element type (1) = linear, etc and number of components
                FE_Q<dim>(heating_degree), 1), // (we have 2 variables: potential and T with 1 component each)
        dof_handler(triangulation), pq(NULL), laplace(NULL), previous_iteration(
        NULL), interp_initial_conditions(false) {
}
// ----------------------------------------------------------------------------------------

template<int dim>
CurrentsAndHeatingStationary<dim>::CurrentsAndHeatingStationary(PhysicalQuantities *pq_, Laplace<dim>* laplace_) :
        ambient_temperature(ambient_temperature_default), fe(FE_Q<dim>(currents_degree), 1, // Finite element type (1) = linear, etc and number of components
                FE_Q<dim>(heating_degree), 1), // (we have 2 variables: potential and T with 1 component each)
        dof_handler(triangulation), pq(pq_), laplace(laplace_), previous_iteration(
        NULL), interp_initial_conditions(false) {
}

template<int dim>
CurrentsAndHeatingStationary<dim>::CurrentsAndHeatingStationary(PhysicalQuantities *pq_, Laplace<dim>* laplace_,
        CurrentsAndHeatingStationary *ch_previous_iteration_) :
        ambient_temperature(ambient_temperature_default), fe(FE_Q<dim>(currents_degree), 1, // Finite element type (1) = linear, etc and number of components
                FE_Q<dim>(heating_degree), 1), // (we have 2 variables: potential and T with 1 component each)
        dof_handler(triangulation), pq(pq_), laplace(laplace_), previous_iteration(
                ch_previous_iteration_), interp_initial_conditions(ch_previous_iteration_ != NULL) {
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::reinitialize(Laplace<dim>* laplace_) {
    reinitialize(laplace_, NULL);
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::reinitialize(Laplace<dim>* laplace_,
        CurrentsAndHeatingStationary *ch_previous_iteration_) {
    dof_handler.clear();
    triangulation.clear();
    interface_map.clear();
    interface_map_field.clear();

    laplace = laplace_;
    previous_iteration = ch_previous_iteration_;

    if (previous_iteration == NULL)
        interp_initial_conditions = false;
    else
        interp_initial_conditions = true;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::set_physical_quantities(PhysicalQuantities *pq_) {
    pq = pq_;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::set_ambient_temperature(double ambient_temperature_) {
    ambient_temperature = ambient_temperature_;
}

template<int dim>
Triangulation<dim>* CurrentsAndHeatingStationary<dim>::get_triangulation() {
    return &triangulation;
}

template<int dim>
DoFHandler<dim>* CurrentsAndHeatingStationary<dim>::get_dof_handler() {
    return &dof_handler;
}

template<int dim>
Vector<double>* CurrentsAndHeatingStationary<dim>::get_solution() {
    return &present_solution;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::import_mesh_from_file(const std::string file_name) {
    MeshPreparer<dim> mesh_preparer;

    mesh_preparer.import_mesh_from_file(&triangulation, file_name);
    mesh_preparer.mark_copper_boundary(&triangulation);

}

template<int dim>
bool CurrentsAndHeatingStationary<dim>::import_mesh_directly(std::vector<Point<dim> > vertices,
        std::vector<CellData<dim> > cells) {
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

    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.mark_copper_boundary(&triangulation);

    return true;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::output_mesh(const std::string file_name) {
    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.output_mesh(&triangulation, file_name);
}

template<int dim>
std::vector<double> CurrentsAndHeatingStationary<dim>::get_temperature(const std::vector<int> &cell_indexes,
        const std::vector<int> &vert_indexes) {

    // Initialise potentials with a value that is immediately visible if it's not changed to proper one
    std::vector<double> temperatures(cell_indexes.size(), 1e15);

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation, 0, cell_indexes[i],
                &dof_handler);

        double temperature = present_solution[dof_cell->vertex_dof_index(vert_indexes[i], 1)];
        temperatures[i] = temperature;
    }
    return temperatures;
}

template<int dim>
std::vector<Tensor<1, dim> > CurrentsAndHeatingStationary<dim>::get_current(
        const std::vector<int> &cell_indexes, const std::vector<int> &vert_indexes) {
    QGauss<dim> quadrature_formula(std::max(currents_degree, heating_degree) + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients);

    std::vector<Tensor<1, dim> > potential_gradients(quadrature_formula.size());
    const FEValuesExtractors::Scalar potential(0);

    std::vector<Tensor<1, dim> > currents(cell_indexes.size());

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation, 0, cell_indexes[i],
                &dof_handler);

        double temperature = present_solution[dof_cell->vertex_dof_index(vert_indexes[i], 1)];

        fe_values.reinit(dof_cell);
        fe_values[potential].get_function_gradients(present_solution, potential_gradients);

        Tensor<1, dim> field = -1.0 * potential_gradients.at(vert_indexes[i]);
        Tensor<1, dim> current = pq->sigma(temperature) * field;

        currents[i] = current;
    }
    return currents;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::setup_system() {
    dof_handler.distribute_dofs(fe);
    //std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
    //      << std::endl;

    DoFRenumbering::component_wise(dof_handler);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    newton_update.reinit(dof_handler.n_dofs());
    present_solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
bool CurrentsAndHeatingStationary<dim>::setup_mapping() {

    double eps = 1e-9;

    // ---------------------------------------------------------------------------------------------
    // Loop over vacuum interface cells

    std::vector<unsigned int> vacuum_interface_indexes;
    std::vector<unsigned int> vacuum_interface_face;
    std::vector<Point<dim> > vacuum_interface_centers;

    typename DoFHandler<dim>::active_cell_iterator vac_cell = laplace->dof_handler.begin_active(),
            vac_endc = laplace->dof_handler.end();
    for (; vac_cell != vac_endc; ++vac_cell) {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
            if (vac_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                vacuum_interface_indexes.push_back(vac_cell->index());
                vacuum_interface_face.push_back(f);
                vacuum_interface_centers.push_back(vac_cell->face(f)->center());
            }
        }
    }
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // Loop over copper interface cells

    typename DoFHandler<dim>::active_cell_iterator cop_cell = dof_handler.begin_active(), cop_endc =
            dof_handler.end();

    for (; cop_cell != cop_endc; ++cop_cell) {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
            if (cop_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                Point<dim> cop_face_center = cop_cell->face(f)->center();
                std::pair<unsigned, unsigned> cop_face_info(cop_cell->index(), f);
                // Loop over vacuum side and find corresponding (cell, face) pair
                for (unsigned int i = 0; i < vacuum_interface_indexes.size(); i++) {
                    if (cop_face_center.distance(vacuum_interface_centers[i]) < eps) {

                        std::pair<unsigned, unsigned> vac_face_info(vacuum_interface_indexes[i],
                                vacuum_interface_face[i]);
                        std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > pair(
                                cop_face_info, vac_face_info);
                        interface_map.insert(pair);
                        break;
                    }
                }
                if (interface_map.count(cop_face_info) == 0)
                    return false;
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    return true;
}

template<int dim>
bool CurrentsAndHeatingStationary<dim>::setup_mapping_field(double smoothing) {

    double eps = 1e-9;

    // ---------------------------------------------------------------------------------------------
    // Loop over vacuum interface cells

    std::vector<Point<dim> > vacuum_interface_centers;
    std::vector<double> vacuum_interface_efield;

    QGauss<dim - 1> face_quadrature_formula(1); // Quadrature with one point
    FEFaceValues<dim> vacuum_fe_face_values(laplace->fe, face_quadrature_formula,
            update_gradients | update_quadrature_points);

    // Electric field values from laplace solver
    std::vector<Tensor<1, dim> > electric_field_value(1);

    typename DoFHandler<dim>::active_cell_iterator vac_cell = laplace->dof_handler.begin_active(),
            vac_endc = laplace->dof_handler.end();
    for (; vac_cell != vac_endc; ++vac_cell) {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
            if (vac_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                // ---
                // Electric field norm in the center (only quadrature point) of the face
                vacuum_fe_face_values.reinit(vac_cell, f);
                vacuum_fe_face_values.get_function_gradients(laplace->solution,
                        electric_field_value);
                double efield_norm = electric_field_value[0].norm();
                // ---

                vacuum_interface_efield.push_back(efield_norm);
                vacuum_interface_centers.push_back(vac_cell->face(f)->center());
            }
        }
    }
    // ---------------------------------------------------------------------------------------------

    // Smoothing: replace a top % with their average + stdev
    if (smoothing > 0.0 && smoothing < 1.0) {
        std::vector<double> efield_vector(vacuum_interface_efield);
        std::sort(efield_vector.begin(), efield_vector.end());
        std::reverse(efield_vector.begin(), efield_vector.end());
        int num_top = smoothing * efield_vector.size();
        std::vector<double> top_elems(efield_vector.begin(), efield_vector.begin() + num_top);
        double mean = vector_mean(top_elems);
        double stdev = vector_stdev(top_elems);
        for (auto &e : vacuum_interface_efield) {
            if (e > mean + stdev) {
                //std::cout << e << " REPLACED BY " << mean + stdev << std::endl;
                e = mean + stdev;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Loop over copper interface cells

    typename DoFHandler<dim>::active_cell_iterator cop_cell = dof_handler.begin_active(), cop_endc =
            dof_handler.end();

    for (; cop_cell != cop_endc; ++cop_cell) {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
            if (cop_cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                Point<dim> cop_face_center = cop_cell->face(f)->center();
                std::pair<unsigned, unsigned> cop_face_info(cop_cell->index(), f);
                // Loop over vacuum side and find corresponding (cell, face) pair
                for (unsigned int i = 0; i < vacuum_interface_centers.size(); i++) {
                    if (cop_face_center.distance(vacuum_interface_centers[i]) < eps) {
                        std::pair<std::pair<unsigned, unsigned>, double> pair(cop_face_info,
                                vacuum_interface_efield[i]);
                        interface_map_field.insert(pair);
                        break;
                    }
                }
                if (interface_map_field.count(cop_face_info) == 0)
                    return false;
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    return true;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::get_surface_nodes(std::vector<Point<dim>>& nodes) {
    const int n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;
    nodes.clear();

    // Loop over copper interface cells
    typename DoFHandler<dim>::active_cell_iterator cell;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (int f = 0; f < n_faces_per_cell; f++)
            if (cell->face(f)->boundary_id() == BoundaryId::copper_surface)
                nodes.push_back(cell->face(f)->center());
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::set_electric_field_bc(const std::vector<double>& elfields) {
    const unsigned n_faces_per_cell = GeometryInfo<dim>::faces_per_cell;

    // Loop over copper interface cells
    typename DoFHandler<dim>::active_cell_iterator cell;
    unsigned i = 0;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        for (unsigned f = 0; f < n_faces_per_cell; ++f)
            if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                std::pair<unsigned, unsigned> face_info(cell->index(), f);
                interface_map_field.insert(
                        std::pair<std::pair<unsigned, unsigned>, double>(face_info, elfields[i++]));
            }
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::set_initial_condition_slow() {
    /* To set the initial condition based on previous solution, we need to find the values
     * at present mesh nodes based on the values of the previous mesh nodes. Present mesh
     * nodes can be outside the old mesh. The number of nodes can also be different. One
     * way to find the values is for every node of the new mesh to loop over all the nodes
     * in the old mesh and find the closest node and set the values from that. This could
     * be sped up by spatially grouping the nodes but this is not done at present moment.
     */

    Triangulation<dim>* previous_triangulation = previous_iteration->get_triangulation();
    DoFHandler<dim>* previous_dof_handler = previous_iteration->get_dof_handler();
    Vector<double>* previous_solution = previous_iteration->get_solution();

    std::vector<Point<dim> > vertex_points;
    std::vector<unsigned> vertex_indexes;

    typename Triangulation<dim>::vertex_iterator it_v = previous_triangulation->begin_vertex(),
            endv = previous_triangulation->end_vertex();
    for (; it_v != endv; it_v++) {
        vertex_points.push_back(it_v->vertex(0));
        vertex_indexes.push_back(it_v->vertex_index(0));
    }
    /* -------------------------------------------------- */

    // Generate a map of "global vertex index" -> "cells" for the old system and the current one
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator> > old_vc_map =
            GridTools::vertex_to_cell_map(*previous_triangulation);

    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator> > vc_map =
            GridTools::vertex_to_cell_map(triangulation);

    /* Iterate all vertices of the current system */
    typename Triangulation<dim>::vertex_iterator it_cs = triangulation.begin_vertex(), end_cs =
            triangulation.end_vertex();
    for (; it_cs != end_cs; it_cs++) {
        const Point<dim> p = it_cs->vertex(0);

        double dist = 1e16;
        unsigned best_i = 0;
        /* Iterate through all vertices of the old system (with the correct x value) */
        /* and find the corresponding vertex index */
        for (unsigned i = 0; i < vertex_points.size(); i++) {
            double c_dist = p.distance(vertex_points[i]);
            if (c_dist < dist) {
                dist = c_dist;
                best_i = i;
            }
        }

        //std::cout << "dist " << dist << std::endl;

        unsigned old_vertex_index = vertex_indexes[best_i];

        typename Triangulation<dim>::active_cell_iterator old_cell =
                *(old_vc_map[old_vertex_index].begin());
        typename DoFHandler<dim>::active_cell_iterator old_dof_cell(previous_triangulation,
                old_cell->level(), old_cell->index(), previous_dof_handler);

        /* Find the potential and temperature values in the corresponding old vertex */
        double pot = 0.0;
        double temp = 0.0;
        for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j) {
            if (old_cell->vertex_index(j) == old_vertex_index) {
                pot = (*previous_solution)[old_dof_cell->vertex_dof_index(j, 0)];
                temp = (*previous_solution)[old_dof_cell->vertex_dof_index(j, 1)];
            }
        }

        typename Triangulation<dim>::active_cell_iterator cell =
                *(vc_map[it_cs->vertex_index(0)].begin());
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation, cell->level(),
                cell->index(), &dof_handler);

        /* Set the potential and temperature values in the new vertex */
        for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j) {
            if (cell->vertex_index(j) == it_cs->vertex_index(0)) {
                present_solution[dof_cell->vertex_dof_index(j, 0)] = pot;
                present_solution[dof_cell->vertex_dof_index(j, 1)] = temp;
            }
        }
    }

}

template<int dim>
void CurrentsAndHeatingStationary<dim>::set_initial_condition() {

    /* If the initial condition is not interpolated from another solution,
     * set temperature at ambient temperature and potential at 0
     */
    if (!interp_initial_conditions) {
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc =
                dof_handler.end();
        for (; cell != endc; ++cell) {
            for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j) {
                present_solution[cell->vertex_dof_index(j, 0)] = 0.0;
                present_solution[cell->vertex_dof_index(j, 1)] = ambient_temperature;
            }
        }
        return;
    }

    /* To set the initial condition based on previous solution, we need to find the values
     * at present mesh nodes based on the values of the previous mesh nodes. Present mesh
     * nodes can be outside the old mesh. The number of nodes can also be different. One
     * way to find the values is for every node of the new mesh to loop over all the nodes
     * in the old mesh and find the closest node and set the values from that. This is too
     * slow. To speed things up the nodes at the previous mesh are grouped near support points.
     */

    Triangulation<dim>* previous_triangulation = previous_iteration->get_triangulation();
    DoFHandler<dim>* previous_dof_handler = previous_iteration->get_dof_handler();
    Vector<double>* previous_solution = previous_iteration->get_solution();

    std::vector<Point<dim> > support_points;

    /* each support point has nearest vertices and their indexes stored */
    std::vector<std::vector<Point<dim> > > sup_to_vertex_point;
    std::vector<std::vector<unsigned> > sup_to_vertex_index;

    unsigned total_vertices = previous_triangulation->n_used_vertices();

    typename Triangulation<dim>::vertex_iterator it_v = previous_triangulation->begin_vertex(),
            endv = previous_triangulation->end_vertex();

    /* Find the bounding box*/
    Point<dim> max_values, min_values;
    for (int d = 0; d < dim; d++) {
        max_values[d] = -1e16;
        min_values[d] = 1e16;
    }
    for (; it_v != endv; it_v++) {
        const Point<dim> p = it_v->vertex(0);
        for (int d = 0; d < dim; d++) {
            if (p[d] > max_values[d])
                max_values[d] = p[d];
            if (p[d] < min_values[d])
                min_values[d] = p[d];
        }
    }

    /* ------------------------------------------------- */
    /* Number of support points in one direction */
    unsigned num_sp_1d = std::pow(std::sqrt(total_vertices), 1. / dim);
    Point<dim> dx;
    for (int d = 0; d < dim; d++)
        dx[d] = (max_values[d] - min_values[d]) / num_sp_1d;

    for (unsigned ix = 0; ix < num_sp_1d; ix++) {
        for (unsigned iy = 0; iy < num_sp_1d; iy++) {
            if (dim == 3) {
                for (unsigned iz = 0; iz < num_sp_1d; iz++)
                    support_points.push_back(
                            Point<dim>(min_values[0] + dx[0] * ix, min_values[1] + dx[1] * iy,
                                    min_values[2] + dx[2] * iz));
            } else {
                support_points.push_back(
                        Point<dim>(min_values[0] + dx[0] * ix, min_values[1] + dx[1] * iy));
            }
        }
    }

    sup_to_vertex_point.resize(support_points.size());
    sup_to_vertex_index.resize(support_points.size());

    it_v = previous_triangulation->begin_vertex();
    for (; it_v != endv; it_v++) {
        Point<dim> p = it_v->vertex(0);

        unsigned index = nearest_point_index(p, support_points);
        sup_to_vertex_point[index].push_back(p);
        sup_to_vertex_index[index].push_back(it_v->vertex_index(0));
    }

    // Remove "empty" support points
    for (unsigned i = 0; i < sup_to_vertex_point.size();) {
        if (sup_to_vertex_point[i].size() == 0) {
            sup_to_vertex_point.erase(sup_to_vertex_point.begin() + i);
            sup_to_vertex_index.erase(sup_to_vertex_index.begin() + i);
            support_points.erase(support_points.begin() + i);
        } else {
            i++;
        }
    }
    // Map of support point -> all nearby vertices and their global indexes is now complete

    /*
     for (int i = 0; i < sup_to_vertex_point.size(); i++) {
     std::cout << sup_to_vertex_point[i].size() << std::endl;
     }
     */
    // Maps from vertex to cells for previous triangulation and the current one
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator> > old_vc_map =
            GridTools::vertex_to_cell_map(*previous_triangulation);
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator> > vc_map =
            GridTools::vertex_to_cell_map(triangulation);

    /* Iterate all vertices of the current system */
    typename Triangulation<dim>::vertex_iterator it_cs = triangulation.begin_vertex(), end_cs =
            triangulation.end_vertex();
    for (; it_cs != end_cs; it_cs++) {
        const Point<dim> p = it_cs->vertex(0);

        // find the nearest support point
        unsigned sp_index = nearest_point_index(p, support_points);

        // find the nearest vertex associated with that support point
        unsigned v_index = nearest_point_index(p, sup_to_vertex_point[sp_index]);

        // global index of the previous system vertex
        unsigned old_vertex_global_index = sup_to_vertex_index[sp_index][v_index];

        typename Triangulation<dim>::active_cell_iterator old_cell =
                *(old_vc_map[old_vertex_global_index].begin());
        typename DoFHandler<dim>::active_cell_iterator old_dof_cell(previous_triangulation,
                old_cell->level(), old_cell->index(), previous_dof_handler);

        /* Find the potential and temperature values in the corresponding old vertex */
        double pot = 0.0;
        double temp = 0.0;
        for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j) {
            if (old_cell->vertex_index(j) == old_vertex_global_index) {
                pot = (*previous_solution)[old_dof_cell->vertex_dof_index(j, 0)];
                temp = (*previous_solution)[old_dof_cell->vertex_dof_index(j, 1)];
            }
        }

        typename Triangulation<dim>::active_cell_iterator cell =
                *(vc_map[it_cs->vertex_index(0)].begin());
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation, cell->level(),
                cell->index(), &dof_handler);

        /* Set the potential and temperature values in the new vertex */
        for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j) {
            if (cell->vertex_index(j) == it_cs->vertex_index(0)) {
                present_solution[dof_cell->vertex_dof_index(j, 0)] = pot;
                present_solution[dof_cell->vertex_dof_index(j, 1)] = temp;
            }
        }
    }
}

// Assembles the linear system for one Newton iteration
template<int dim>
void CurrentsAndHeatingStationary<dim>::assemble_system_newton() {

    TimerOutput timer(std::cout, TimerOutput::never, TimerOutput::wall_times);
    timer.enter_section("Pre-assembly");

    QGauss<dim> quadrature_formula(std::max(currents_degree, heating_degree) + 1);
    QGauss<dim - 1> face_quadrature_formula(
            std::max(std::max(currents_degree, heating_degree), laplace->shape_degree) + 1);

    FEValues<dim> fe_values(fe, quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
            update_values | update_gradients | update_normal_vectors | update_quadrature_points
                    | update_JxW_values);

    // For evaluating E field at the copper surface
    FEFaceValues<dim> vacuum_fe_face_values(laplace->fe, face_quadrature_formula,
            update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // ---------------------------------------------------------------------------------------------
    // Declaring variables once at the start (for performance reasons)

    // The previous solution values in the cell quadrature points
    std::vector<Tensor<1, dim>> prev_sol_potential_gradients(n_q_points);
    std::vector<double> prev_sol_temperature_values(n_q_points);
    std::vector<Tensor<1, dim>> prev_sol_temperature_gradients(n_q_points);

    // The previous solution values in the face quadrature points
    std::vector<Tensor<1, dim>> prev_sol_face_potential_gradients(n_face_q_points);
    std::vector<double> prev_sol_face_temperature_values(n_face_q_points);
    std::vector<Tensor<1, dim>> prev_sol_face_temperature_gradients(n_face_q_points);

    // Shape function values and gradients (arrays for every cell DOF)
    std::vector<double> potential_phi(dofs_per_cell);
    std::vector<Tensor<1, dim> > potential_phi_grad(dofs_per_cell);
    std::vector<double> temperature_phi(dofs_per_cell);
    std::vector<Tensor<1, dim> > temperature_phi_grad(dofs_per_cell);

    // Electric field values from laplace solver
    std::vector<Tensor<1, dim> > electric_field_values(n_face_q_points);
    // ----------------------------------------------------------------------------------------------

    const FEValuesExtractors::Scalar potential(0);
    const FEValuesExtractors::Scalar temperature(1);

    timer.exit_section();
    timer.enter_section("Loop header");

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc =
            dof_handler.end();
    for (; cell != endc; ++cell) {

        fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs = 0;

        fe_values[potential].get_function_gradients(present_solution, prev_sol_potential_gradients);
        fe_values[temperature].get_function_values(present_solution, prev_sol_temperature_values);
        fe_values[temperature].get_function_gradients(present_solution,
                prev_sol_temperature_gradients);

        timer.exit_section();
        timer.enter_section("Matrix assembly 1");

        // ---------------------------------------------------------------------------------------------
        // Local matrix assembly
        // ---------------------------------------------------------------------------------------------
        for (unsigned int q = 0; q < n_q_points; ++q) {

            double prev_temp = prev_sol_temperature_values[q];
            const Tensor<1, dim> prev_pot_grad = prev_sol_potential_gradients[q];
            const Tensor<1, dim> prev_temp_grad = prev_sol_temperature_gradients[q];

            double sigma = pq->sigma(prev_temp);
            double dsigma = pq->dsigma(prev_temp);
            double kappa = pq->kappa(prev_temp);
            double dkappa = pq->dkappa(prev_temp);

            for (unsigned int k = 0; k < dofs_per_cell; ++k) {
                potential_phi_grad[k] = fe_values[potential].gradient(k, q);
                temperature_phi[k] = fe_values[temperature].value(k, q);
                temperature_phi_grad[k] = fe_values[temperature].gradient(k, q);
            }

            timer.exit_section();
            timer.enter_section("Matrix assembly 2");

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    cell_matrix(i, j) += (-(potential_phi_grad[i] * sigma * potential_phi_grad[j])
                            - (potential_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
                            + (temperature_phi[i] * 2 * sigma * prev_pot_grad
                                    * potential_phi_grad[j])
                            + (temperature_phi[i] * dsigma * prev_pot_grad * prev_pot_grad
                                    * temperature_phi[j])
                            - (temperature_phi_grad[i] * dkappa * prev_temp_grad
                                    * temperature_phi[j])
                            - (temperature_phi_grad[i] * kappa * temperature_phi_grad[j]))
                            * fe_values.JxW(q);
                    cell_matrix(i, j) += (-(potential_phi_grad[i] * sigma * potential_phi_grad[j])
                            - (potential_phi_grad[i] * dsigma * prev_pot_grad * temperature_phi[j])
                            + (temperature_phi[i] * 2 * sigma * prev_pot_grad
                                    * potential_phi_grad[j])
                            + (temperature_phi[i] * dsigma * prev_pot_grad * prev_pot_grad
                                    * temperature_phi[j])
                            - (temperature_phi_grad[i] * dkappa * prev_temp_grad
                                    * temperature_phi[j])
                            - (temperature_phi_grad[i] * kappa * temperature_phi_grad[j]))
                            * fe_values.JxW(q);
                }
                cell_rhs(i) += (potential_phi_grad[i] * sigma * prev_pot_grad
                        - temperature_phi[i] * sigma * prev_pot_grad * prev_pot_grad
                        + temperature_phi_grad[i] * kappa * prev_temp_grad) * fe_values.JxW(q);
            }
            timer.exit_section();
            timer.enter_section("Matrix assembly 1");
        }
        // ---------------------------------------------------------------------------------------------
        timer.exit_section();
        timer.enter_section("Rhs assembly");
        // ---------------------------------------------------------------------------------------------
        // Local right-hand side assembly
        // ---------------------------------------------------------------------------------------------
        // integration over the boundary (cell faces)
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {

                if (cell->face(f)->boundary_id() == BoundaryId::copper_surface) {
                    fe_face_values.reinit(cell, f);

                    fe_face_values[potential].get_function_gradients(present_solution,
                            prev_sol_face_potential_gradients);
                    fe_face_values[temperature].get_function_values(present_solution,
                            prev_sol_face_temperature_values);
                    fe_face_values[temperature].get_function_gradients(present_solution,
                            prev_sol_face_temperature_gradients);

                    // ---------------------------------------------------------------------------------------------
                    // Vacuum side stuff
                    // find the corresponding vacuum side face to the copper side face
                    std::pair<unsigned, unsigned> cop_cell_info = std::pair<unsigned, unsigned>(
                            cell->index(), f);
                    // check if the corresponding vacuum face exists in our mapping
                    assert(interface_map_field.count(cop_cell_info) == 1);
                    /*
                     std::pair<unsigned, unsigned> vac_cell_info = interface_map[cop_cell_info];
                     // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
                     typename DoFHandler<dim>::active_cell_iterator vac_cell(&(laplace->triangulation),
                     0, vac_cell_info.first, &(laplace->dof_handler));

                     vacuum_fe_face_values.reinit(vac_cell, vac_cell_info.second);
                     vacuum_fe_face_values.get_function_gradients(laplace->solution, electric_field_values);
                     */
                    double e_field = interface_map_field[cop_cell_info];
                    // ---------------------------------------------------------------------------------------------

                    // loop through the quadrature points
                    for (unsigned int q = 0; q < n_face_q_points; ++q) {

                        double prev_temp = prev_sol_face_temperature_values[q];
                        const Tensor<1, dim> prev_pot_grad = prev_sol_face_potential_gradients[q];
                        const Tensor<1, dim> prev_temp_grad = prev_sol_face_temperature_gradients[q];

                        const Tensor<1, dim> normal_vector = fe_face_values.normal_vector(q);

                        double dsigma = pq->dsigma(prev_temp);
                        double dkappa = pq->dkappa(prev_temp);
                        //double e_field = electric_field_values[q].norm();
                        double emission_current = pq->emission_current(e_field, prev_temp);
                        // Nottingham heat flux in
                        // (eV*A/nm^2) -> (eV*n*q_e/(s*nm^2)) -> (J*n/(s*nm^2)) -> (W/nm^2)
                        double nottingham_flux = -1.0 * pq->nottingham_de(e_field, prev_temp)
                                * emission_current;

                        for (unsigned int k = 0; k < dofs_per_cell; ++k) {
                            potential_phi[k] = fe_face_values[potential].value(k, q);
                            temperature_phi[k] = fe_face_values[temperature].value(k, q);
                        }
                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            cell_rhs(i) += (-(potential_phi[i] * emission_current)
                                    - (temperature_phi[i] * nottingham_flux))
                                    * fe_face_values.JxW(q);

                            /*
                             std::cout << (potential_phi[i] * sigma * normal_vector * prev_pot_grad) << " "
                             << (potential_phi[i] * emission_current) << std::endl;
                             */
                            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                                cell_matrix(i, j) += ((potential_phi[i] * normal_vector * dsigma
                                        * prev_pot_grad * temperature_phi[j])
                                        + (temperature_phi[i] * normal_vector * dkappa
                                                * prev_temp_grad * temperature_phi[j]))
                                        * fe_face_values.JxW(q);
                            }
                        }
                    }
                }
            }
        }
        // ---------------------------------------------------------------------------------------------

        timer.exit_section();
        timer.enter_section("Loop footer");

        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
        timer.exit_section();
        timer.enter_section("Loop header");
    }

    timer.exit_section();
    timer.enter_section("Post assembly");

    // Setting Dirichlet boundary values //

    // 0 potential at the bulk bottom boundary
    std::map<types::global_dof_index, double> current_dirichlet;
    VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
            ZeroFunction<dim>(2), current_dirichlet, fe.component_mask(potential));

    MatrixTools::apply_boundary_values(current_dirichlet, system_matrix, newton_update, system_rhs);

    // Set 0 temperature BC, as the initial condition already has correct dirichlet BCs
    std::map<types::global_dof_index, double> temperature_dirichlet;
    VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
            ZeroFunction<dim>(2), temperature_dirichlet, fe.component_mask(temperature));

    /*
     if (first_iteration) {
     VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
     ConstantFunction<dim>(ambient_temperature, 2), temperature_dirichlet, fe.component_mask(temperature));
     } else {
     VectorTools::interpolate_boundary_values(dof_handler, BoundaryId::copper_bottom,
     ZeroFunction<dim>(2), temperature_dirichlet, fe.component_mask(temperature));
     }
     */

    MatrixTools::apply_boundary_values(temperature_dirichlet, system_matrix, newton_update,
            system_rhs);

    timer.exit_section();
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::solve() {
    // CG doesn't work as the matrix is not symmetric

    // GMRES
    /*
     SolverControl solver_control(400000, 1e-9);
     SolverGMRES<> solver_gmres(solver_control, SolverGMRES<>::AdditionalData(50));
     solver_gmres.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
     std::cout << "   " << solver_control.last_step() << " GMRES iterations needed to obtain convergence." << std::endl;
     */

    // UMFPACK solver
    deallog << "Solving linear system with UMFPACK... " << std::endl;
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(newton_update, system_rhs);
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::run() {

    std::cout << "/---------------------------------------------------------------/" << std::endl
            << "CurrentsAndHeating run():" << std::endl;

    double temperature_tolerance = 1.0;

    Timer timer;

    setup_system();
    setup_mapping_field();

    // Sets the initial state
    // Dirichlet BCs need to hold for this state
    // and 0 dirichlet BC should be applied for all Newton iterations
    set_initial_condition();

    std::cout << "    Setup and IC: " << timer.wall_time() << " s" << std::endl;

    // Newton iterations
    for (unsigned int iteration = 0; iteration < 5; ++iteration) {
        std::cout << "/--------------------------------/" << std::endl;
        std::cout << "Newton iteration " << iteration << std::endl;

        timer.restart();
        // reset the state of the linear system
        system_matrix.reinit(sparsity_pattern);
        system_rhs.reinit(dof_handler.n_dofs());
        std::cout << "    Reset state: " << timer.wall_time() << " s" << std::endl;
        timer.restart();

        timer.restart();
        assemble_system_newton();

        std::cout << "    Assembly: " << timer.wall_time() << " s" << std::endl;
        timer.restart();

        solve();
        present_solution.add(2.0, newton_update); // alpha = 1.0

        std::cout << "    Solver: " << timer.wall_time() << " s" << std::endl;
        timer.restart();

        output_results("output/solution.vtk", iteration);
        std::cout << "    output_results: " << timer.wall_time() << " s" << std::endl;
        timer.restart();

        std::cout << "    ||u_k-u_{k-1}||_L2 = " << newton_update.l2_norm() << std::endl;
        std::cout << "    ||u_k-u_{k-1}||_Linf = " << newton_update.linfty_norm() << std::endl;
        std::cout << "    Residual = " << system_rhs.l2_norm() << std::endl;

        if (newton_update.linfty_norm() < temperature_tolerance) {
            std::cout << "    Maximum temperature change less than tolerance: converged!"
                    << std::endl;
            break;
        }
    }
    std::cout << "/---------------------------------------------------------------/" << std::endl;
}

template<int dim>
double CurrentsAndHeatingStationary<dim>::run_specific(double temperature_tolerance, int max_newton_iter,
        bool file_output, std::string out_fname, bool print, double sor_alpha,
        double ic_interp_treshold, bool skip_field_mapping) {

    if (pq == NULL || laplace == NULL
            || (interp_initial_conditions && previous_iteration == NULL)) {
        std::cerr << "Error: pointer uninitialized! Exiting temperature calculation..."
                << std::endl;
        return -1.0;
    }
    if ((laplace->solution).size() == 0) {
        std::cerr << "Error: Laplace solution hasn't been calculation." << std::endl;
        return -1.0;
    }

    if (interp_initial_conditions) {
        double prev_max_temp = (previous_iteration->present_solution).linfty_norm();
        if (prev_max_temp > ic_interp_treshold) {
            if (print)
                std::cout << "        Interpolating initial conditions" << std::endl;
        } else {
            interp_initial_conditions = false;
            if (print)
                std::cout
                        << "        Using default initial conditions (peak temp. lower than threshold)"
                        << std::endl;
        }
    } else if (print)
        std::cout << "        Using default initial conditions" << std::endl;

    Timer timer;

    if (!skip_field_mapping) {
        if (!setup_mapping_field()) {
            std::cerr
                    << "Error: Couldn't make a correct mapping between copper and vacuum faces on the interface."
                    << "Make sure that the face elements have one-to-one correspondence there."
                    << std::endl;
            return -1.0;
        }
    }

    if (print) {
        printf("        Mapping setup done, time: %.2f\n", timer.wall_time());
        timer.restart();
    }

    // Sets the initial state
    // Dirichlet BCs need to hold for this state
    // and 0 dirichlet BC should be applied for all Newton iterations
    set_initial_condition();

    if (print) {
        printf("        Initial condition setup, time: %.2f\n", timer.wall_time());
        timer.restart();
    }

    // Output initial condition
    if (file_output) {
        output_results(out_fname, 0);
        if (print) {
            printf("        Output initial condition, time: %.2f\n", timer.wall_time());
            timer.restart();
        }
    }

    double temperature_error = 1e15;

    // Newton iterations
    for (int iteration = 1; iteration < max_newton_iter + 1; ++iteration) {

        system_matrix.reinit(sparsity_pattern);
        system_rhs.reinit(dof_handler.n_dofs());

        // Set dirichlet BSs as 0, as they're already set in  the initial condition
        assemble_system_newton();
        double assemble_time = timer.wall_time();
        timer.restart();

        solve();
        present_solution.add(sor_alpha, newton_update);
        double solution_time = timer.wall_time();
        timer.restart();

        double max_temp = present_solution.linfty_norm();
        if (max_temp > temperature_stopping_condition) {
            std::cerr
                    << "WARNING: Peak temperature surpassed the stopping condition "
                            + to_string(temperature_stopping_condition)
                    << "... Stopping calculation." << std::endl;
            break;
        }

        if (file_output)
            output_results(out_fname, iteration);
        double output_time = timer.wall_time();
        timer.restart();

        temperature_error = newton_update.linfty_norm();

        double potential_rel_error = 0.0;
        for (unsigned int i = 0; i<newton_update.size()/2; i++) {
            double local_rel_error = std::abs(sor_alpha*newton_update[i]/present_solution[i]);
            if (local_rel_error > potential_rel_error) potential_rel_error = local_rel_error;
        }

        if (print) {
            printf("        iter: %2d; t_error: %7.3f; p_rel_err: %2.0e; assemble_time: %5.2f;"
                   " sol_time: %5.2f; outp_time: %5.2f\n",
                   iteration, temperature_error, potential_rel_error, assemble_time, solution_time, output_time);
        }

        if (temperature_error < temperature_tolerance && potential_rel_error < 0.5) {
            break;
        }
    }
    return temperature_error;
}

template<int dim>
void CurrentsAndHeatingStationary<dim>::output_results(const std::string file_name,
        const int iteration) const {
    std::string file_name_mod = file_name;

    if (iteration >= 0) {
        const int start = file_name.find_last_of('.');
        const int end = file_name.size();
        std::string ext = ".vtk";
        if (start != -1) {
            ext = file_name.substr(start, end);
            file_name_mod = file_name.substr(0, start);
        }
        file_name_mod += "-" + std::to_string(iteration) + ext;
    }
    std::vector<std::string> solution_names { "potential", "temperature" };

    CurrentPostProcessorStat<dim> current_post_processor(pq); // needs to be before data_out
    SigmaPostProcessorStat<dim> sigma_post_processor(pq); // needs to be before data_out

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(present_solution, solution_names);
    data_out.add_data_vector(present_solution, current_post_processor);
    data_out.add_data_vector(present_solution, sigma_post_processor);
    data_out.build_patches();

    try {
        std::ofstream output(file_name_mod.c_str());
        data_out.write_vtk(output);
    } catch (...) {
        std::cerr << "WARNING: Couldn't open " + file_name_mod << ". ";
        std::cerr << "Output is not saved." << std::endl;
    }
}

template class CurrentsAndHeatingStationary<2> ;
template class CurrentsAndHeatingStationary<3> ;

} // namespace fch

