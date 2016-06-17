/*
 * DealII.cpp
 /*
 *  Created on: 11.2.2016
 *      Author: veske
 */

#include "DealII.h"
#include <fstream>

using namespace std;
using namespace dealii;

namespace femocs {

// Laplace solver constructor
// Number in fe() determines the interpolation type. 1 is linear etc.
DealII::DealII() :
        fe(POLY_DEGREE), dof_handler(triangulation), efield(0), neumann(0) {
}

const void DealII::set_neumann(const double efield) {
    // The field on upper boundary is on the limit of zmax->inf exactly 27x higher than efield
    this->neumann = efield / 27.0;
    this->efield = -1.0*efield;
}

// Get long range electric field value
const double DealII::get_efield() {
    return efield;
}

// Get number of degrees of freedom
const int DealII::get_n_dofs() {
    return dof_handler.n_dofs();
}

// Get number of used nodes in mesh
const int DealII::get_n_nodes() {
    return triangulation.n_used_vertices();
}

// Generate simple mesh for test purposes
const void DealII::make_simple_mesh() {
    Triangulation<DIM> tr1, tr2;

    vector<Point<DIM> > vertices1(DIM + 1);
    vector<Point<DIM> > vertices2(DIM + 1);
    vector<Point<DIM> > vertices3(DIM + 1);

    vertices1[0] = Point<DIM>(1.0, 0.0, 0.7);
    vertices1[1] = Point<DIM>(-1.0, 0.0, 0.7);
    vertices1[2] = Point<DIM>(0.0, 1.0, -0.7);
    vertices1[3] = Point<DIM>(0.0, -1.0, -0.7);

    vertices2[0] = Point<DIM>(1.0 + 2.0, 0.0 + 2.0, 0.7 + 2.0);
    vertices2[1] = Point<DIM>(-1.0 + 2.0, 0.0 + 2.0, 0.7 + 2.0);
    vertices2[2] = Point<DIM>(0.0 + 2.0, 1.0 + 2.0, -0.7 + 2.0);
    vertices2[3] = Point<DIM>(0.0 + 2.0, -1.0 + 2.0, -0.7 + 2.0);

    vertices3[0] = Point<DIM>(1.0 + 4.0, 0.0 + 4.0, 0.7 + 4.0);
    vertices3[1] = Point<DIM>(-1.0 + 4.0, 0.0 + 4.0, 0.7 + 4.0);
    vertices3[2] = Point<DIM>(0.0 + 4.0, 1.0 + 4.0, -0.7 + 4.0);
    vertices3[3] = Point<DIM>(0.0 + 4.0, -1.0 + 4.0, -0.7 + 4.0);

    GridGenerator::simplex(tr1, vertices1);
    GridGenerator::simplex(tr2, vertices2);
    GridGenerator::merge_triangulations(tr1, tr2, triangulation);

    tr1.clear();
    tr2.clear();

    GridGenerator::simplex(tr1, vertices3);
    GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
}

// Import mesh from vtk or msh file
const void DealII::import_file(const string file_name) {
    const string file_type = get_file_type(file_name);
    expect(file_type == "vtk" || file_type == "msh", "Unknown file type!")

    GridIn<DIM, DIM> gi;
    gi.attach_triangulation(triangulation);
    ifstream in_file(file_name);

    if (file_type == "vtk")
        gi.read_vtk(in_file);
    else if (file_type == "msh")
        gi.read_msh(in_file);
}

// Modified version of grid_in function,
// https://github.com/dealii/dealii/blob/master/source/grid/grid_in.cc
const void DealII::import_tethex_mesh(tethex::Mesh* mesh) {
    const unsigned int n_verts = mesh->get_n_vertices();
    const unsigned int n_elems = mesh->get_n_hexahedra();
    const unsigned int n_faces = mesh->get_n_quadrangles();
    unsigned int i;

    vector<Point<DIM> > vertices(n_verts);       // array for vertices
    vector<CellData<DIM> > cells(n_elems);       // array for elements
    SubCellData subcelldata;
    subcelldata.boundary_quads.reserve(n_faces); // array for faces

    // copy vertices
    for (unsigned int vertex = 0; vertex < n_verts; ++vertex)
        for (i = 0; i < DIM; ++i)
            vertices[vertex](i) = mesh->get_vertex(vertex).get_coord(i);

    // copy quadrangles
    for (unsigned int face = 0; face < n_faces; ++face) {
        subcelldata.boundary_quads[face] = CellData<2>();
        for (i = 0; i < n_verts_per_face; ++i)
            subcelldata.boundary_quads[face].vertices[i] = mesh->get_quadrangle(face).get_vertex(i);
    }

    // copy hexahedra
    for (unsigned int elem = 0; elem < n_elems; ++elem) {
        cells[elem] = CellData<DIM>();
        for (i = 0; i < n_verts_per_elem; ++i)
            cells[elem].vertices[i] = mesh->get_hexahedron(elem).get_vertex(i);
    }

    // Check consistency of subcelldata
    Assert(subcelldata.check_consistency(DIM), ExcInternalError());
    // Do some clean-up on vertices...
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    // ... and on cells
    GridReordering<DIM, DIM>::invert_all_cells_of_negative_grid(vertices, cells);
//    GridReordering<DIM, DIM>::reorder_cells(cells);

    triangulation.create_triangulation_compatibility(vertices, cells, subcelldata);
}

// Import tetrahedral mesh
const void DealII::import_tetgen_mesh(femocs::Mesh* mesh) {
    Triangulation<DIM> tr1, tr2;
    vector<Point<DIM>> vertices1(DIM + 1), vertices2(DIM + 1);
    SimpleElement selem;
    Point3 node;

    int i, j, elem;
    int n_elems = mesh->get_n_elems();
    triangulation.clear();

    if (n_elems < 1) return;

    selem = mesh->get_simpleelem(0);
    for (i = 0; i < n_verts_per_elem; ++i) {
        node = mesh->get_node(selem[i]);
        vertices1[i] = Point<DIM>(node.x, node.y, node.z);
    }

    if (n_elems == 1) {
        GridGenerator::simplex(triangulation, vertices1);
        return;
    }

    selem = mesh->get_simpleelem(1);
    for (i = 0; i < n_verts_per_elem; ++i) {
        node = mesh->get_node(selem[i]);
        vertices2[i] = Point<DIM>(node.x, node.y, node.z);
    }

    GridGenerator::simplex(tr1, vertices1);
    GridGenerator::simplex(tr2, vertices2);
    GridGenerator::merge_triangulations(tr1, tr2, triangulation);
    tr1.clear();
    tr2.clear();

    // loop through tetrahedra, convert them into hexahedra and add to big triangulations
    for (elem = 2; elem < n_elems; ++elem) {
        if (elem % 10 == 0) cout << elem << "/" << n_elems << "\n";
        selem = mesh->get_simpleelem(elem);
        for (i = 0; i < n_verts_per_elem; ++i) {
            node = mesh->get_node(selem[i]);
            vertices1[i] = Point<DIM>(node.x, node.y, node.z);
        }

        GridGenerator::simplex(tr1, vertices1);
        GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
        tr1.clear();
    }
}

// Mark the boundary faces of mesh
const void DealII::mark_boundary(const AtomReader::Sizes* sizes, const AtomReader::Types* types) {
    typename Triangulation<DIM>::active_face_iterator face;
    const double eps = 0.1;

    // Loop through the faces and mark them according the location of their centre
    for (face = triangulation.begin_face(); face != triangulation.end_face(); ++face)
        if (face->at_boundary()) {
            if (on_boundary(face->center()[0], sizes->xmin, sizes->xmax, eps))
                face->set_all_boundary_ids(types->type_edge);
            else if (on_boundary(face->center()[1], sizes->ymin, sizes->ymax, eps))
                face->set_all_boundary_ids(types->type_edge);
            else if (on_boundary(face->center()[2], sizes->zmaxbox, eps))
                face->set_all_boundary_ids(types->type_zmax);
            else
                face->set_all_boundary_ids(types->type_surf);
        }
}

// Setup initial grid and number the vertices i.e. distribute degrees of freedom.
const void DealII::setup_system() {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());

    laplace_solution.reinit(dof_handler.n_dofs());
}

// Insert boundary conditions to the system
const void DealII::assemble_system(const AtomReader::Types* types) {
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
    const unsigned int n_dofs_per_elem = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    // Declare cell matrix and right-hand-side matrix (sets coordinate system so we get real positive cells)
    FullMatrix<double> cell_matrix(n_dofs_per_elem, n_dofs_per_elem);
    Vector<double> cell_rhs(n_dofs_per_elem);

    // Create a vector of local degrees of freedom for each cell.
    vector<types::global_dof_index> local_dof_indices(n_dofs_per_elem);

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
            for (i = 0; i < n_dofs_per_elem; ++i)
                cell_rhs(i) += fe_values.shape_value(i, q_index) * space_charge
                        * fe_values.JxW(q_index);

            // Assemble the matrix
            for (i = 0; i < n_dofs_per_elem; ++i)
                for (j = 0; j < n_dofs_per_elem; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q_index)
                            * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index);

            // Cycle for faces of each cell
            for (unsigned int f = 0; f < n_faces_per_elem; ++f)
                // Check if face is located at top boundary
                if (cell->face(f)->boundary_id() == types->type_zmax) {
                    fe_face_values.reinit(cell, f);
                    // Set boundary conditions
                    for (unsigned int fq_index = 0; fq_index < n_face_q_points; ++fq_index)
                        for (i = 0; i < n_dofs_per_elem; ++i)
                            // Set Neumann boundary value
                            cell_rhs(i) += fe_face_values.shape_value(i, fq_index) * neumann
                                    * fe_face_values.JxW(fq_index);
                }
        }

        // Apply set conditions and rewrite system matrix and rhs
        cell->get_dof_indices(local_dof_indices);
        for (i = 0; i < n_dofs_per_elem; ++i)
            for (j = 0; j < n_dofs_per_elem; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

        for (i = 0; i < n_dofs_per_elem; ++i)
            system_rhs(local_dof_indices[i]) += cell_rhs(i);

    }

    // Declare boundaries
    map<types::global_dof_index, double> copper_boundary_value;

    // Add Dirichlet' boundary condition to faces denoted as surface
    VectorTools::interpolate_boundary_values(dof_handler, types->type_surf, ZeroFunction<DIM>(),
            copper_boundary_value);

    // Apply boundary values to system matrix
    MatrixTools::apply_boundary_values(copper_boundary_value, system_matrix, laplace_solution,
            system_rhs);
}

// Run the calculation with conjugate gradient solver
const void DealII::solve_cg() {
    // Declare Conjugate Gradient solver tolerance and max number of iterations
    SolverControl solver_control(10000, 1e-9);
    SolverCG<> solver(solver_control);
    solver.solve(system_matrix, laplace_solution, system_rhs, PreconditionIdentity());
}

// Run the calculation with UMFPACK solver
const void DealII::solve_umfpack() {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(laplace_solution, system_rhs);
}

// Calculate electric field from arbitrary point inside the mesh
const Tensor<1, DIM> DealII::get_elfield_at_point(Point<DIM> &point) {
    return -1.0 * VectorTools::point_gradient(dof_handler, laplace_solution, point);
}

// Calculate electric field at node point
const Tensor<1, DIM> DealII::get_elfield_at_node(const int &cell_indx, const int &vert_indx) {
    vector<Tensor<1, DIM>> solution_gradients;
    QTrapez<DIM> only_vertices_quadrature_formula;
    FEValues<DIM> fe_values(fe, only_vertices_quadrature_formula, update_gradients);
    solution_gradients.resize(only_vertices_quadrature_formula.size());

    typename DoFHandler<DIM>::active_cell_iterator cell;

    cell = dof_handler.begin_active();
    for (unsigned int i = 0; i < cell_indx; ++i)
        ++cell;

    fe_values.reinit(cell);
    fe_values.get_function_gradients(laplace_solution, solution_gradients);

    return -1.0 * solution_gradients.at(vert_indx);
}

// Calculate potential at arbitrary point inside the mesh
const double DealII::get_potential_at_point(Point<DIM> &point) {
    return VectorTools::point_value(dof_handler, laplace_solution, point);
}

// Get potential at mesh node
const double DealII::get_potential_at_node(const int &cell_indx, const int &vert_indx) {
    typename DoFHandler<DIM>::active_cell_iterator cell;

    cell = dof_handler.begin_active();
    for (unsigned int i = 0; i < cell_indx; ++i)
        ++cell;

    return laplace_solution(cell->vertex_dof_index(vert_indx, 0));
}

// Write the potential and electric field to the file
const void DealII::output_results(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string ftype = get_file_type(file_name);
    expect(ftype == "vtk" || ftype == "eps", "Unsupported file type!");

    LaplacePostProcessor field_calculator("Electric_field");
    DataOut<DIM> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(laplace_solution, "Potential");
    data_out.add_data_vector(laplace_solution, field_calculator);
    data_out.build_patches();

    ofstream output(file_name);
    if (ftype == "vtk")
        data_out.write_vtk(output);
    else if (ftype == "eps") data_out.write_eps(output);
}

// Write the mesh to file
const void DealII::output_mesh(const string file_name) {
#if not DEBUGMODE
    return;
#endif
    string ftype = get_file_type(file_name);
    expect(ftype == "vtk" || ftype == "msh" || ftype == "eps", "Unsupported file type!");

    ofstream outfile(file_name);
    GridOut gout;
    gout.write_msh(triangulation, outfile);
    if (ftype == "vtk")
        gout.write_vtk(triangulation, outfile);
    else if (ftype == "msh")
        gout.write_msh(triangulation, outfile);
    else if (ftype == "eps") gout.write_eps(triangulation, outfile);
}

} /* namespace femocs */
