/*
 * mesh_input.h
 *
 *  Created on: Jul 27, 2016
 *      Author: kristjan
 */

#ifndef INCLUDE_MESH_PREPARER_H_
#define INCLUDE_MESH_PREPARER_H_

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

namespace fch {

using namespace std;
using namespace dealii;

enum BoundaryId {
    copper_surface = 1,
    vacuum_top = 2,
    copper_bottom = 3,
    vacuum_sides = 4,
    copper_sides = 5
};

/**
 * Sets boundary identifications for the mesh boundaries
 */
template<int dim>
class MeshPreparer {

    const std::string get_file_ext(const std::string file_name);

public:

    /**
     * Import hexahedral (3D) or quadrangular (2D) mesh from file
     * @param triangulation  pointer to mesh where the data is imported to
     * @param file_name name of the input mesh file (.vtk or .msh format)
     */
    void import_mesh_from_file(Triangulation<dim> *triangulation,
            const std::string file_name);

    /**
     * Writes mesh to vtk file
     * @param triangulation  pointer to any mesh data
     * @param name           name of the file
     */
    void output_mesh(Triangulation<dim> *triangulation, std::string name);

    /**
     * Sets boundary identification to vacuum mesh top and lateral sides and to the copper-vacuum boundary
     * @param triangulation  pointer to vacuum mesh data
     */
    void mark_vacuum_boundary(Triangulation<dim> *triangulation);

    /**
      * Sets boundary identification to bulk mesh bottom and lateral sides and to the copper-vacuum boundary
      * @param triangulation  pointer to bulk mesh data
      */
    void mark_copper_boundary(Triangulation<dim> *triangulation);

    /**
     * @TODO figure out description
     */
    void mark_top_and_bottom_boundary(Triangulation<dim> *triangulation);

    /**
     * Delete cells with given id from the mesh
     * @param triangulation  pointer to the mesh data
     * @param id  id of the cells to be deleted
     * @return  pointer to the cleaned mesh data
     */
    Triangulation<dim> remove_cells_with_id(Triangulation<dim> *triangulation, int id);

    /**
     * Marks boundary (top, bottom, sides, vacuum-copper surface) faces (3D) or edges (2D) in the mesh
     * @param triangulation  pointer to the mesh data
     * @param top     id for top faces|edges
     * @param bottom  id for bottom faces|edges
     * @param sides   id for lateral faces|edges
     * @param other   id for vacuum-copper surface faces|edges
     */
    void mark_boundary(Triangulation<dim> *triangulation, char top, char bottom, char sides,
            char other);

};

} // namespace fch

#endif /* INCLUDE_MESH_PREPARER_H_ */
