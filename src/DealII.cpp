/*
 * DealII.cpp
 /*
 *  Created on: 11.2.2016
 *      Author: veske
 */

#include "DealII.h"
#include <fstream>

#include <omp.h>

using namespace std;
using namespace dealii;

namespace femocs {

// Laplace solver constructor
// Number in fe() determines the interpolation type. 1 is linear etc.
DealII::DealII() :
        fe(POLY_DEGREE), dof_handler(triangulation), neumann(0) {
}

// Specify the Neumann boundary condition value
const void DealII::set_neumann(const double neumann) {
    this->neumann = neumann;
}

// Get long range electric field value
const double DealII::get_elfield() {
    return -1.0 * neumann;
}

// Get number of degrees of freedom
const int DealII::get_n_dofs() {
    return dof_handler.n_dofs();
}

// Get number of used nodes in mesh
const int DealII::get_n_nodes() {
    return triangulation.n_used_vertices();
}

// Get number of used faces in mesh
const int DealII::get_n_faces() {
    return triangulation.n_active_faces();
}

// Get number of used cells in mesh
const int DealII::get_n_cells() {
    return triangulation.n_active_cells();
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
    expect(file_type == "vtk" || file_type == "msh", "Unknown file type: " + file_type)

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

//    triangulation.create_triangulation_compatibility(vertices, cells, subcelldata);
    triangulation.create_triangulation_compatibility(vertices, cells, SubCellData());
}

const void DealII::import_tethex_mesh_vol2(tethex::Mesh* mesh) {
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

    // copy hexahedra
    for (unsigned int elem = 0; elem < n_elems; ++elem) {
        cells[elem] = CellData<DIM>();
        for (i = 0; i < n_verts_per_elem; ++i)
            cells[elem].vertices[i] = mesh->get_hexahedron(elem).get_vertex(i);
    }

    // Do some clean-up on vertices...
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    // ... and on cells
    GridReordering<DIM, DIM>::invert_all_cells_of_negative_grid(vertices, cells);
//    GridReordering<DIM, DIM>::reorder_cells(cells);

    triangulation.create_triangulation_compatibility(vertices, cells, SubCellData());
}

// Import tetrahedral mesh
// WARNING: It's extremely slow function
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
        // show progress after every 10th step
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
const void DealII::mark_boundary(const AtomReader::Sizes* sizes) {


    typename Triangulation<DIM>::active_face_iterator face;
    const double eps = 0.1;

    // Loop through the faces and mark them according the location of their centre
    for (face = triangulation.begin_face(); face != triangulation.end_face(); ++face)
        if (face->at_boundary()) {
            if (on_boundary(face->center()[0], sizes->xmin, sizes->xmax, eps))
                face->set_all_boundary_ids(TYPES.EDGE);
            else if (on_boundary(face->center()[1], sizes->ymin, sizes->ymax, eps))
                face->set_all_boundary_ids(TYPES.EDGE);
            else if (on_boundary(face->center()[2], sizes->zmaxbox, eps))
                face->set_all_boundary_ids(TYPES.ZMAX);
            else
                face->set_all_boundary_ids(TYPES.SURFACE);
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
const void DealII::assemble_system() {
    // Set up quadrature system for quads and faces
    const QGauss<DIM> quadrature_formula(2);
    const QGauss<DIM-1> face_quadrature_formula(2);
    unsigned int i, j, q_index;

    // Calculate necessary values (derived from weak form of Laplace equation)
    FEValues<DIM> fe_values(fe, quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<DIM> fe_face_values(fe, face_quadrature_formula,
            update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Parametrize necessary entities.
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
//    const double space_charge = 0.0;

    // Declare cell matrix and right-hand-side matrix (sets coordinate system so we get real positive cells)
    FullMatrix<double> cell_matrix(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double> cell_rhs(n_dofs_per_cell);

    // Create a vector of local degrees of freedom for each cell.
    vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

    //Start iterator cycles over all cells to set local terms for each cell
    DoFHandler<DIM>::active_cell_iterator cell;
    for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        // Integration of the cell by looping over all quadrature points
        for (q_index = 0; q_index < quadrature_formula.size(); ++q_index) {
//            // Assemble the right hand side by integrating over the shape function i times the
//            // right hand side function; we choose it to be the function with constant value one
//            for (i = 0; i < n_dofs_per_cell; ++i)
//                cell_rhs(i) += fe_values.shape_value(i, q_index) * space_charge
//                    * fe_values.JxW(q_index);

            // Assemble the matrix
            for (i = 0; i < n_dofs_per_cell; ++i)
                for (j = 0; j < n_dofs_per_cell; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q_index) * fe_values.JxW(q_index)
                        * fe_values.shape_grad(j, q_index);
        }

        // Cycle for faces of each cell
        for (unsigned int f = 0; f < n_faces_per_elem; ++f) {
            // Check if face is located at top boundary
            if (cell->face(f)->boundary_id() == TYPES.ZMAX) {
                fe_face_values.reinit(cell, f);

                // Apply Neumann boundary condition
                for (q_index = 0; q_index < face_quadrature_formula.size(); ++q_index)
                    for (i = 0; i < n_dofs_per_cell; ++i)
                        cell_rhs(i) += fe_face_values.shape_value(i, q_index) * neumann
                                * fe_face_values.JxW(q_index);
            }
        }

        // Apply set conditions and rewrite system matrix & rhs
        cell->get_dof_indices(local_dof_indices);
        for (i = 0; i < n_dofs_per_cell; ++i)
            for (j = 0; j < n_dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

        for (i = 0; i < n_dofs_per_cell; ++i)
            system_rhs(local_dof_indices[i]) += cell_rhs(i);

    }

    // Declare boundaries
    map<types::global_dof_index, double> bv;

    // Add Dirichlet' boundary condition to faces denoted as surface
    VectorTools::interpolate_boundary_values(dof_handler, TYPES.SURFACE, ZeroFunction<DIM>(), bv);

    // Apply boundary values to system matrix
    MatrixTools::apply_boundary_values(bv, system_matrix, laplace_solution, system_rhs);
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

// Calculate electric field at arbitrary point inside the mesh
const Tensor<1, DIM> DealII::get_elfield_at_point(Point<DIM> &point) {
    return -1.0 * VectorTools::point_gradient(dof_handler, laplace_solution, point);
}

// Calculate electric field at arbitrary point inside the mesh
const Tensor<1, DIM> DealII::get_elfield_at_point(const double x, const double y, const double z) {
    return -1.0 * VectorTools::point_gradient(dof_handler, laplace_solution, Point<DIM>(x, y, z));
}

// Calculate electric field at a mesh node
const Tensor<1, DIM> DealII::get_elfield_at_node(const int &cell_indx, const int &vert_indx) {
    vector<Tensor<1, DIM>> solution_gradients;
    QTrapez<DIM> only_vertices_quadrature_formula;
    FEValues<DIM> fe_values(fe, only_vertices_quadrature_formula, update_gradients);

    solution_gradients.resize(only_vertices_quadrature_formula.size());

    typename DoFHandler<DIM>::active_cell_iterator cell;
    cell = dof_handler.begin_active();

    std::advance(cell, cell_indx);
    fe_values.reinit(cell);
    fe_values.get_function_gradients(laplace_solution, solution_gradients);

    return -1.0 * solution_gradients.at(vert_indx);
}

// Calculate electric field at a set of mesh nodes
const vector<Vec3> DealII::get_elfield_at_node(const vector<int> &cell_indxs, const vector<int> &vert_indxs) {
    const int n_cells = cell_indxs.size();
    vector<Vec3> elfield(n_cells);

    vector<Tensor<1, DIM>> solution_gradients;
    QTrapez<DIM> only_vertices_quadrature_formula;
    FEValues<DIM> fe_values(fe, only_vertices_quadrature_formula, update_gradients);
    solution_gradients.resize(only_vertices_quadrature_formula.size());

    // Generate vector with indices [0, n_atoms-1]
    vector<int> sort_indxs(n_cells);
    size_t n(0);
    generate(sort_indxs.begin(), sort_indxs.end(), [&]{ return n++; });

    // Sort indexes by the cell_indxs to make it possible to iterate sequentially though the cells
    auto comparator = [&cell_indxs](int a, int b){ return cell_indxs[a] < cell_indxs[b]; };
    sort(sort_indxs.begin(), sort_indxs.end(), comparator);

    int i, j;
    int si = sort_indxs[0];
    typename DoFHandler<DIM>::active_cell_iterator cell;

    // Iterate through all the cells and get the electric field from ones listed in cell_indxs
    for (i = 0, j = 0, cell = dof_handler.begin_active(); cell != dof_handler.end(); i++, cell++)
        if (i == cell_indxs[si]) {
            fe_values.reinit(cell);
            fe_values.get_function_gradients(laplace_solution, solution_gradients);

            Tensor<1, DIM> ef = -1.0 * solution_gradients.at(vert_indxs[si]);
            elfield[si] = Vec3(ef[0], ef[1], ef[2]);

            si = sort_indxs[++j];
        }

    return elfield;
}

// Calculate potential at arbitrary point inside the mesh
const double DealII::get_potential_at_point(Point<DIM> &point) {
    return VectorTools::point_value(dof_handler, laplace_solution, point);
}

// Calculate potential at arbitrary point inside the mesh
const double DealII::get_potential_at_point(const double x, const double y, const double z) {
    return VectorTools::point_value(dof_handler, laplace_solution, Point<DIM>(x, y, z));
}

// Get potential at a mesh node
const double DealII::get_potential_at_node(const int &cell_indx, const int &vert_indx) {
    typename DoFHandler<DIM>::active_cell_iterator cell;
    cell = dof_handler.begin_active();
    std::advance(cell, cell_indx);

    return laplace_solution(cell->vertex_dof_index(vert_indx, 0));
}

// Get potential at a set of mesh nodes
const vector<double> DealII::get_potential_at_node(const vector<int> &cell_indxs, const vector<int> &vert_indxs) {
    const int n_cells = cell_indxs.size();

    vector<double> potentials(n_cells);

    // Generate vector with indices [0, n_atoms-1]
    vector<int> sort_indxs(n_cells);
    size_t n(0);
    generate(sort_indxs.begin(), sort_indxs.end(), [&]{ return n++; });

    // Sort indexes by the cell_indxs to make it possible to iterate sequentially though the cells
    auto comparator = [&cell_indxs](int a, int b){ return cell_indxs[a] < cell_indxs[b]; };
    sort(sort_indxs.begin(), sort_indxs.end(), comparator);

    int i, j;
    int si = sort_indxs[0];
    typename DoFHandler<DIM>::active_cell_iterator cell;

    // Iterate through all the cells and get the potential from the ones listed in cell_indxs
    for (i = 0, j = 0, cell = dof_handler.begin_active(); cell != dof_handler.end(); i++, cell++)
        if (i == cell_indxs[si]) {
            potentials[si] = laplace_solution( cell->vertex_dof_index(vert_indxs[si], 0) );
            si = sort_indxs[++j];
        }

    return potentials;
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
