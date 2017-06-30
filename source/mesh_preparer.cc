/*
 * mesh_preparer.cc
 *
 *  Created on: Jul 27, 2016
 *      Author: kristjan
 */

#include "mesh_preparer.h"

#include <fstream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace fch {

template<int dim>
const string MeshPreparer<dim>::get_file_ext(const string file_name) {
    const int start = file_name.find_last_of('.') + 1;
    const int end = file_name.size();
    if (start > end)
        return "";
    return file_name.substr(start, end);
}

template<int dim>
void MeshPreparer<dim>::output_mesh(Triangulation<dim> *triangulation, std::string file_name) {
    const string file_type = get_file_ext(file_name);
    std::ofstream out(file_name);
    GridOut grid_out;
    grid_out.set_flags(GridOutFlags::Msh(true, true));

    try {
        if (file_type == "msh") {
            grid_out.write_msh(*triangulation, out);
        } else {
            grid_out.write_vtk(*triangulation, out);
        }
    } catch (...) {
        std::cerr << "WARNING: Couldn't open " + file_name << ". ";
        std::cerr << "Mesh is not saved." << std::endl;
    }
}

template<int dim>
void MeshPreparer<dim>::import_mesh_from_file(Triangulation<dim> *triangulation,
        const string file_name) {
    const string file_type = get_file_ext(file_name);
    if (!(file_type == "vtk" || file_type == "msh")) {
        cout << "Error: Unknown file type!" << endl;
        return;
    }

    GridIn<dim, dim> gi;

    gi.attach_triangulation(*triangulation);
    ifstream in_file(file_name);

    if (file_type == "vtk")
        gi.read_vtk(in_file);
    else if (file_type == "msh")
        gi.read_msh(in_file);
}

template<int dim>
void MeshPreparer<dim>::mark_boundary(Triangulation<dim> *triangulation, char top, char bottom,
        char sides, char other) {
    typename Triangulation<dim>::active_face_iterator face;
    double eps = 0.1;
    double xmax = -1e16, ymax = -1e16, zmax = -1e16;
    double xmin = 1e16, ymin = 1e16, zmin = 1e16;

    // Loop through the faces and find maximum and minimum values for coordinates
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
        if (face->at_boundary()) {
            double x = face->center()[0];
            double y = face->center()[1];
            if (x > xmax)
                xmax = x;
            if (x < xmin)
                xmin = x;
            if (y > ymax)
                ymax = y;
            if (y < ymin)
                ymin = y;

            if (dim == 3) {
                double z = face->center()[2];
                if (z > zmax)
                    zmax = z;
                if (z < zmin)
                    zmin = z;
            }
        }
    }
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
        if (face->at_boundary()) {
            if (dim == 2) {
                double x = face->center()[0];
                double y = face->center()[1];

                if (fabs(x - xmin) < eps || fabs(x - xmax) < eps) {
                    face->set_all_boundary_ids(sides);
                } else if (fabs(y - ymax) < eps) {
                    face->set_all_boundary_ids(top);
                } else if (fabs(y - ymin) < eps) {
                    face->set_all_boundary_ids(bottom);
                } else {
                    face->set_all_boundary_ids(other);
                }
            } else if (dim == 3) {
                double x = face->center()[0];
                double y = face->center()[1];
                double z = face->center()[2];

                if (fabs(x - xmin) < eps || fabs(x - xmax) < eps || fabs(y - ymin) < eps
                        || fabs(y - ymax) < eps) {
                    face->set_all_boundary_ids(sides);
                } else if (fabs(z - zmax) < eps) {
                    face->set_all_boundary_ids(top);
                } else if (fabs(z - zmin) < eps) {
                    face->set_all_boundary_ids(bottom);
                } else {
                    face->set_all_boundary_ids(other);
                }
            }

        }
    }
}

template<int dim>
void MeshPreparer<dim>::mark_top_and_bottom_boundary(Triangulation<dim> *triangulation) {
    mark_boundary(triangulation, BoundaryId::vacuum_top, BoundaryId::copper_bottom,
            BoundaryId::vacuum_sides, BoundaryId::vacuum_sides);
}

template<int dim>
void MeshPreparer<dim>::mark_vacuum_boundary(Triangulation<dim> *triangulation) {
    mark_boundary(triangulation, BoundaryId::vacuum_top, BoundaryId::copper_surface,
            BoundaryId::vacuum_sides, BoundaryId::copper_surface);
}

template<int dim>
void MeshPreparer<dim>::mark_copper_boundary(Triangulation<dim> *triangulation) {
    mark_boundary(triangulation, BoundaryId::copper_surface, BoundaryId::copper_bottom,
            BoundaryId::copper_sides, BoundaryId::copper_surface);
}

template<int dim>
Triangulation<dim> MeshPreparer<dim>::remove_cells_with_id(Triangulation<dim> *triangulation,
        int id) {
    std::set<typename Triangulation<dim>::active_cell_iterator> removal;
    typename Triangulation<dim>::active_cell_iterator cell;
    for (cell = triangulation->begin_active(); cell != triangulation->end(); ++cell) {
        if (cell->material_id() == id)
            removal.insert(cell);
    }
    Triangulation<dim> tmp_triang;
    GridGenerator::create_triangulation_with_removed_cells(*triangulation, removal, tmp_triang);

    return tmp_triang;
}

// Give all template possibilities here (alternative is to write implementation in .h file)
template class MeshPreparer<2> ;
template class MeshPreparer<3> ;

} // namespace fch
