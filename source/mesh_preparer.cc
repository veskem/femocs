/*
 * mesh_preparer.cc
 *
 *  Created on: Jul 27, 2016
 *      Author: kristjan
 */


#include "mesh_preparer.h"

#include <fstream>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

template <int dim>
void MeshPreparer<dim>::output_mesh(Triangulation<dim> *triangulation, std::string name) {
	std::ofstream out(name);
	GridOut grid_out;
	grid_out.write_vtk(*triangulation, out);
	std::cout << "Grid written to " << name << std::endl;
}

template <int dim>
const string MeshPreparer<dim>::get_file_ext(const string file_name) {
    const int start = file_name.find_last_of('.') + 1;
    const int end = file_name.size();
    return file_name.substr(start, end);
}

template <int dim>
void MeshPreparer<dim>::import_mesh_from_file(Triangulation<dim> *triangulation, const string file_name) {
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

template <int dim>
void MeshPreparer<dim>::mark_vacuum_boundary(Triangulation<dim> *triangulation) {
	typename Triangulation<dim>::active_face_iterator face;
	double eps = 0.1;
	double xmax = -1e16, ymax = -1e16, zmax = -1e16;
	double xmin = 1e16, ymin = 1e16;

    // Loop through the faces and find maximum and minimum values for coordinates
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
    	if (face->at_boundary()) {
			double x = face->center()[0];
			double y = face->center()[1];
			double z = face->center()[2];
			if (x > xmax) xmax = x;
			if (x < xmin) xmin = x;
			if (y > ymax) ymax = y;
			if (y < ymin) ymin = y;
			if (z > zmax) zmax = z;
    	}
    }
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
    	if (face->at_boundary()) {
    		double x = face->center()[0];
			double y = face->center()[1];
			double z = face->center()[2];

			if (fabs(x-xmin) < eps || fabs(x-xmax) < eps || fabs(y-ymin) < eps || fabs(y-ymax) < eps) {
				face->set_all_boundary_ids(BoundaryId::vacuum_sides);
			} else if (fabs(z-zmax) < eps) {
				face->set_all_boundary_ids(BoundaryId::vacuum_top);
			} else {
				face->set_all_boundary_ids(BoundaryId::copper_surface);
			}
    	}
    }
}

template <int dim>
void MeshPreparer<dim>::mark_copper_boundary(Triangulation<dim> *triangulation) {
	typename Triangulation<dim>::active_face_iterator face;
	double eps = 0.1;
	double xmax = -1e16, ymax = -1e16;
	double xmin = 1e16, ymin = 1e16, zmin = 1e16;

    // Loop through the faces and find maximum and minimum values for coordinates
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
    	if (face->at_boundary()) {
			double x = face->center()[0];
			double y = face->center()[1];
			double z = face->center()[2];
			if (x > xmax) xmax = x;
			if (x < xmin) xmin = x;
			if (y > ymax) ymax = y;
			if (y < ymin) ymin = y;
			if (z < zmin) zmin = z;
    	}
    }
    for (face = triangulation->begin_face(); face != triangulation->end_face(); ++face) {
    	if (face->at_boundary()) {
    		double x = face->center()[0];
			double y = face->center()[1];
			double z = face->center()[2];

			if (fabs(x-xmin) < eps || fabs(x-xmax) < eps || fabs(y-ymin) < eps || fabs(y-ymax) < eps) {
				face->set_all_boundary_ids(BoundaryId::copper_sides);
			} else if (fabs(z-zmin) < eps) {
				face->set_all_boundary_ids(BoundaryId::copper_bottom);
			} else {
				face->set_all_boundary_ids(BoundaryId::copper_surface);
			}
    	}
    }
}

// Give all template possibilities here (alternative is to write implementation in .h file)
template class MeshPreparer<2>;
template class MeshPreparer<3>;

/*
template void MeshPreparer::draw_mesh<2>(Triangulation<2>&, std::string);
template void MeshPreparer::import_mesh_from_file<2>(Triangulation<2> *triangulation, const string file_name);

template void MeshPreparer::draw_mesh<3>(Triangulation<3>&, std::string);
template void MeshPreparer::import_mesh_from_file<3>(Triangulation<3> *triangulation, const string file_name);
*/
