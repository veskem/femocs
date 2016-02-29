/*
 * DealII.cpp
 *
 *  Created on: 11.2.2016
 *      Author: veske
 */

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "DealII.h"

using namespace std;
using namespace dealii;
namespace femocs {

// Define main class constructor. Number in fe(2) determines the interpolation type. 1 would be linear etc.
//DealII::DealII() : fe(1), dof_handler(triangulation) {
DealII::DealII(const int poly_degree, const double neumann, Femocs::SimuCell* simucell) :
        fe(poly_degree), dof_handler(triangulation) {
    this->neumann = neumann;
    types.bulk = simucell->type_bulk;
    types.side = simucell->type_edge;
    types.surface = simucell->type_surf;
    types.top = simucell->type_zmax;
    types.vacuum = simucell->type_vacuum;
    sizes.xmax = simucell->xmax;
    sizes.xmin = simucell->xmin;
    sizes.ymax = simucell->ymax;
    sizes.ymin = simucell->ymin;
    sizes.zmax = simucell->zmax;
    sizes.zmin = simucell->zmin;
    sizes.zmaxbox = simucell->zmaxbox;
    sizes.zminbox = simucell->zminbox;
}

// Import mesh from vtk or msh file
void DealII::import_file(const string file_name) {
    GridIn<DIM, DIM> gi;

    gi.attach_triangulation(triangulation);
    string file_type = get_file_type(file_name);
    ifstream in_file(file_name);

    if (file_type == "vtk")
        gi.read_vtk(in_file);
    else if (file_type == "msh")
        gi.read_msh(in_file);
    else
        cout << "Unknown file type!\n";
}

// Extract the file type from file name
const string DealII::get_file_type(const string& file_name) {
    int start = file_name.find_last_of('.') + 1;
    int end = file_name.size();
    return file_name.substr(start, end);
}

// Import hexahedral mesh from tethex
void DealII::import_tethex_mesh(tethex::Mesh* mesh) {
    int n_elems, n_faces, n_verts;
    int n_verts_in_face, n_verts_in_elem;

    n_verts = mesh->get_n_vertices();
    n_faces = mesh->get_n_quadrangles();
    n_elems = mesh->get_n_hexahedra();
    n_verts_in_face = GeometryInfo<DIM-1>::vertices_per_cell;
    n_verts_in_elem = GeometryInfo<DIM>::vertices_per_cell;

    if (n_elems < 0) {
        cout << "Number of elements in input mesh < 1!";
        return;
    }

/*
    // magic numbers that turn elements to follow right sequence
    const int order_of_faces[] = { 0, 1, 3, 2, 4, 5, 7, 6 };
    // generate hexahedra vertice indexes
    vector<CellData<2>> faces(n_faces, CellData<2>());

    for (int face = 0; face < n_faces; ++face) {
        faces[face].material_id = 0;
        for (int vert = 0; vert < n_verts_in_face; ++vert)
            faces[face].vertices[vert] = face * n_verts_in_face + order_of_faces[vert];
    }
 //*/   
    // magic numbers that turn elements to follow right sequence
    const int order_of_elems[] = { 0, 1, 3, 2, 4, 5, 7, 6 };
    // generate hexahedra vertice indexes
    vector<CellData<DIM>> elems(n_elems, CellData<DIM>());

    for (int elem = 0; elem < n_elems; ++elem) {
        elems[elem].material_id = 0;
        for (int vert = 0; vert < n_verts_in_elem; ++vert)
            elems[elem].vertices[vert] = elem * n_verts_in_elem + order_of_elems[vert];
    }

    // magic numbers that turn vertices to follow right sequence
    const int order_of_verts[] = { 0, 1, 5, 4, 3, 2, 6, 7 };
    vector<Point<DIM>> vertices(n_elems * n_verts_in_elem);
    int vert_cntr = 0;

    // transfer nodes from tethex to deal.ii mesh
    for (int elem = 0; elem < n_elems; ++elem)
        for (int vert = 0; vert < n_verts_in_elem; ++vert) {
            int vert2 = mesh->get_hexahedron(elem).get_vertex(order_of_verts[vert]);
            double v1 = mesh->get_vertex(vert2).get_coord(0);
            double v2 = mesh->get_vertex(vert2).get_coord(1);
            double v3 = mesh->get_vertex(vert2).get_coord(2);
            vertices[vert_cntr++] = Point<DIM>(v1, v2, v3);
        }

    //   invert_all_cells_of_negative_grid and reorder_cells function of GridReordering
    //   before creating the triangulation
    //   const bool  use_new_style_ordering = true;
    //   GridReordering<DIM>::invert_all_cells_of_negative_grid (vertices, elems);
    //   GridReordering<DIM>::reorder_cells(elems, use_new_style_ordering);

    // combine vertices and hexahedra into triangulation object
//    triangulation.create_triangulation(vertices, faces, SubCellData());
    triangulation.create_triangulation(vertices, elems, SubCellData());

//    Triangulation<2> tri_2d;
//    tri_2d.create_triangulation(vertices, faces, SubCellData());
//
//    GridGenerator::merge_triangulations(tri_2d, triangulation, triangulation);
//
//   tri_2d.clear();
}

void DealII::import_tethex_mesh_old(tethex::Mesh* mesh) {
    int n_elems, n_verts, n_verts_in_elem;

    n_elems = mesh->get_n_hexahedra();

    if (n_elems < 0) {
        cout << "Number of elements in input mesh < 1!";
        return;
    }

    n_verts = mesh->get_n_vertices();
    n_verts_in_elem = GeometryInfo<DIM>::vertices_per_cell; //mesh->get_hexahedron(0).get_n_vertices();

    // copy nodes from tethex to deal.ii mesh
    vector<Point<DIM>> vertices(n_verts);
    for (int vert = 0; vert < n_verts; ++vert) {
        double v1 = mesh->get_vertex(vert).get_coord(0);
        double v2 = mesh->get_vertex(vert).get_coord(1);
        double v3 = mesh->get_vertex(vert).get_coord(2);
        vertices[vert] = Point<DIM>(v1, v2, v3);
    }

    // copy hexahedra from tethex to deal.ii mesh
    vector<CellData<DIM>> elems(n_elems, CellData<DIM>());
    for (int elem = 0; elem < n_elems; ++elem) {
        elems[elem].material_id = 0;
        for (int vert = 0; vert < n_verts_in_elem; ++vert)
            elems[elem].vertices[vert] = mesh->get_hexahedron(elem).get_vertex(vert);
    }

    // invert_all_cells_of_negative_grid and reorder_cells function of GridReordering 
    // before creating the triangulation

    const bool use_new_style_ordering = true;
    GridReordering<DIM>::invert_all_cells_of_negative_grid(vertices, elems);
    GridReordering<DIM>::reorder_cells(elems, use_new_style_ordering);

    triangulation.create_triangulation(vertices, elems, SubCellData());
}

// Import tetrahedral mesh
void DealII::import_tetgen_mesh(shared_ptr<Mesh> mesh) {
    const int n_nodes_in_elem = 4;
    const int n_coords = 3;

    Triangulation<DIM> tr1;
    Triangulation<DIM> tr2;
    vector<Point<DIM>> vertices1(3 + 1);
    vector<Point<DIM>> vertices2(3 + 1);

    int i, j, node, elem;
    int n_elems = mesh->getNelems();
    triangulation.clear();

    if (n_elems < 1) return;

    for (i = 0; i < n_nodes_in_elem; ++i) {
        node = mesh->getElem(0, i);
        for (j = 0; j < n_coords; ++j)
            vertices1[i](j) = mesh->getNode(node, j);
    }

    if (n_elems == 1) {
        GridGenerator::simplex(triangulation, vertices1);
        return;
    }

    for (i = 0; i < n_nodes_in_elem; ++i) {
        node = mesh->getElem(1, i);
        for (j = 0; j < n_coords; ++j)
            vertices2[i](j) = mesh->getNode(node, j);
    }
    GridGenerator::simplex(tr1, vertices1);
    GridGenerator::simplex(tr2, vertices2);
    GridGenerator::merge_triangulations(tr1, tr2, triangulation);
    tr1.clear();
    tr2.clear();

    // loop through tetrahedra, convert them into hexahedron and add to big triangulations
    for (elem = 2; elem < n_elems; ++elem) {
        if (elem % 10 == 0) cout << elem << "/" << n_elems << "\n";
        for (i = 0; i < n_nodes_in_elem; ++i) {
            node = mesh->getElem(elem, i);
            for (j = 0; j < n_coords; ++j)
                vertices1[i](j) = mesh->getNode(node, j);
        }
        GridGenerator::simplex(tr1, vertices1);
        GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
        tr1.clear();
    }
}

// Generate simple mesh for test purposes
void DealII::make_simple_mesh() {
    const unsigned int d = 3;
    Triangulation<d> tr1;
    Triangulation<d> tr2;

    vector<Point<d> > vertices1(d + 1);
    vector<Point<d> > vertices2(d + 1);
    vector<Point<d> > vertices3(d + 1);
    vector<Point<d> > vertices4(d + 1);

    vertices1[0](0) = 1.;
    vertices1[0](1) = 0.;
    vertices1[0](2) = .7;
    vertices1[1](0) = -1.;
    vertices1[1](1) = 0.;
    vertices1[1](2) = .7;
    vertices1[2](0) = 0.;
    vertices1[2](1) = 1.;
    vertices1[2](2) = -.7;
    vertices1[3](0) = 0.;
    vertices1[3](1) = -1.;
    vertices1[3](2) = -.7;

    vertices2[0](0) = 1. + 2.0;
    vertices2[0](1) = 0. + 2.0;
    vertices2[0](2) = .7 + 2.0;
    vertices2[1](0) = -1. + 2.0;
    vertices2[1](1) = 0. + 2.0;
    vertices2[1](2) = .7 + 2.0;
    vertices2[2](0) = 0. + 2.0;
    vertices2[2](1) = 1. + 2.0;
    vertices2[2](2) = -.7 + 2.0;
    vertices2[3](0) = 0. + 2.0;
    vertices2[3](1) = -1. + 2.0;
    vertices2[3](2) = -.7 + 2.0;

    vertices3[0](0) = 1. + 4.0;
    vertices3[0](1) = 0. + 4.0;
    vertices3[0](2) = .7 + 4.0;
    vertices3[1](0) = -1. + 4.0;
    vertices3[1](1) = 0. + 4.0;
    vertices3[1](2) = .7 + 4.0;
    vertices3[2](0) = 0. + 4.0;
    vertices3[2](1) = 1. + 4.0;
    vertices3[2](2) = -.7 + 4.0;
    vertices3[3](0) = 0. + 4.0;
    vertices3[3](1) = -1. + 4.0;
    vertices3[3](2) = -.7 + 4.0;

    vertices4[0](0) = 1. + 1.0;
    vertices4[0](1) = 0. + 1.0;
    vertices4[0](2) = .7 + 1.0;
    vertices4[1](0) = -1. + 0.0;
    vertices4[1](1) = 0. + 0.0;
    vertices4[1](2) = .7 + 0.0;
    vertices4[2](0) = 0. + 0.0;
    vertices4[2](1) = 1. + 0.0;
    vertices4[2](2) = -.7 + 0.0;
    vertices4[3](0) = 0. + 0.0;
    vertices4[3](1) = -1. + 0.0;
    vertices4[3](2) = -.7 + 0.0;

    GridGenerator::simplex(tr1, vertices1);
    GridGenerator::simplex(tr2, vertices2);
    GridGenerator::merge_triangulations(tr1, tr2, triangulation);

    tr1.clear();
    tr2.clear();

    GridGenerator::simplex(tr1, vertices3);
    GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
    tr1.clear();

    GridGenerator::simplex(tr1, vertices4);
    GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
}

// Setup initial grid and number the vertices i.e. distribute degrees of freedom.
void DealII::setup_system() {
    dof_handler.distribute_dofs(fe);

    cout << "\nNumber of degrees of freedom: " << dof_handler.n_dofs() << endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());

    solution.reinit(dof_handler.n_dofs());
}

void DealII::distort_solution(const double dist_ampl) {
    for (int i = 0; i < solution.size(); ++i)
        solution[i] += dist_ampl * random();
}

void DealII::distort_solution_const(const double dist_ampl) {
    for (int i = 0; i < solution.size(); ++i)
        solution[i] += dist_ampl;
}

void DealII::distort_solution_one(const double dist_ampl, const int i) {
    solution[i] += dist_ampl;
}

// Mark the boundary faces of mesh
void DealII::mark_boundary() {
/*
    typename Triangulation<DIM>::active_cell_iterator cell;
    const unsigned int faces_per_cell = GeometryInfo<DIM>::faces_per_cell;

    // Loop through the faces and mark them according the location of its centre
    for(cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
        for (unsigned int f = 0; f < faces_per_cell; ++f)
            if(cell->face(f)->at_boundary()) {
                if (on_boundary(cell->face(f)->center()[0], sizes.xmin, sizes.xmax))
                    cell->face(f)->set_all_boundary_ids(types.side);
                else if (on_boundary(cell->face(f)->center()[1], sizes.ymin, sizes.ymax))
                    cell->face(f)->set_all_boundary_ids(types.side);
                else if (on_boundary(cell->face(f)->center()[2], sizes.zmaxbox, sizes.zmaxbox))
                    cell->face(f)->set_all_boundary_ids(types.top);
                else
                    cell->face(f)->set_all_boundary_ids(types.surface);
            }
//*/
//*
    typename Triangulation<DIM>::active_face_iterator face;

   
    // Loop through the faces and mark them according the location of its centre
    for(face = triangulation.begin_active_face(); face != triangulation.end(); ++face) {
        if(face->at_boundary()) {
            if (on_boundary(face->center()[0], sizes.xmin, sizes.xmax))
                face->set_all_boundary_ids(types.side);
            else if (on_boundary(face->center()[1], sizes.ymin, sizes.ymax))
                face->set_all_boundary_ids(types.side);
            else if (on_boundary(face->center()[2], sizes.zmaxbox, sizes.zmaxbox))
                face->set_all_boundary_ids(types.top);
            else
                face->set_all_boundary_ids(types.surface);
        }
    }

 //*/
}

// Function to determine whether the center of face is on the boundary of simulation cell or not
bool DealII::on_boundary(const double face, const double face_min, const double face_max) {
    const double eps = 0.1;
    bool b1 = fabs(face - face_min) < eps;
    bool b2 = fabs(face_max - face) < eps;
    return b1 || b2;
}

// Insert boundary conditions to the system
void DealII::assemble_system() {
    // Set up quadrature system for quads and faces
    const QGauss<3> quadrature_formula(DIM);
    const QGauss<2> face_quadrature_formula(DIM);
    unsigned int i, j;


    // Calculate necessary values (derived from weak form of Laplace equation)
    FEValues<DIM> fe_values(fe, quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<DIM> fe_face_values(fe, face_quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Parametrize necessary entities.
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int faces_per_cell = GeometryInfo<DIM>::faces_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    // Declare cell matrix and right-hand-side matrix (sets coordinate system so we get real positive cells)
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    // Create a vector of local degrees of freedom for each cell.
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    //Start iterator cycles over all cells to set local terms for each cell
    DoFHandler<DIM>::active_cell_iterator cell;

    const double space_charge = 0.0;

    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        // Integration of the cell by looping over all quadrature points
        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) {
            // Assemble the right hand side by integrating over the shape function i times the
            // right hand side function; we choose it to be the function with constant value one
            for (i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += fe_values.shape_value(i, q_index) * space_charge * fe_values.JxW(q_index);

            // Assemble the matrix
            for (i = 0; i < dofs_per_cell; ++i)
                for (j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q_index)
                            * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index);

            // Cycle for faces of each cell
            for (unsigned int f = 0; f < faces_per_cell; ++f)
                // Check if face is located at top boundary  
                if ( cell->face(f)->boundary_id() == types.top ) {
                    fe_face_values.reinit(cell, f);
                    // Set boundary conditions
                    for (unsigned int fq_index = 0; fq_index < n_face_q_points; ++fq_index)
                        for (i = 0; i < dofs_per_cell; ++i)
                            // Set Neumann boundary value
                            cell_rhs(i) += fe_face_values.shape_value(i, fq_index) * neumann
                                    * fe_face_values.JxW(fq_index);
                }
        }

        // Apply set conditions and rewrite system matrix and rhs
        cell->get_dof_indices(local_dof_indices);
        for (i = 0; i < dofs_per_cell; ++i)
            for (j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

        for (i = 0; i < dofs_per_cell; ++i)
            system_rhs(local_dof_indices[i]) += cell_rhs(i);

    }

    // Declare boundaries
    map<types::global_dof_index, double> copper_boundary_value;

    // Add Dirichlet' boundary condition to faces denoted as surface
    VectorTools::interpolate_boundary_values(dof_handler, types.surface, ZeroFunction<DIM>(), copper_boundary_value);

    // Apply boundary values to system matrix
    MatrixTools::apply_boundary_values(copper_boundary_value, system_matrix, solution, system_rhs);
}

// Run the calculation with conjugate gradient solver
void DealII::solve_cg() {
    //declare Conjugate Gradient solver tolerance and number if iterations used to achieve a converged solution
    SolverControl solver_control(10000, 1e-9);
    SolverCG<> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    /*
     // Initialize solver
     PreconditionSSOR<> preconditioner;
     preconditioner.initialize(system_matrix, 1.2);

     //run solver
     solver.solve(system_matrix, solution, system_rhs, preconditioner);

     //make constraints apply to solution
     constraints.distribute(solution);
     //*/
}

// Run the calculation with UMFPACK solver
void DealII::solve_umfpack() {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
}

// Output the calculation results to vtk file
void DealII::output_results(const string file_name) {
    DataOut<DIM> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output(file_name);
    data_out.write_vtk(output);
}

// Output the mesh 
void DealII::output_mesh(const string file_name) {
    ofstream outfile(file_name);
    GridOut gout;
    gout.write_vtk(triangulation, outfile);
}

} /* namespace femocs */
