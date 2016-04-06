/*
 * DealII.cpp
 /*
 *  Created on: 11.2.2016
 *      Author: veske
 */

#include "DealII.h"

using namespace std;
using namespace dealii;

namespace femocs {

// Define main class constructor. Number in fe(2) determines the interpolation type. 1 would be linear etc.
//DealII::DealII() : fe(1), dof_handler(triangulation) {
DealII::DealII(const int poly_degree, const double neumann) :
        fe(poly_degree), dof_handler(triangulation) {
    this->neumann = neumann;
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

// Inspiration from the grid_in source,
// https://github.com/dealii/dealii/blob/master/source/grid/grid_in.cc
void DealII::import_tethex_mesh(tethex::Mesh* mesh) {

    const unsigned int n_verts_per_elem = GeometryInfo<3>::vertices_per_cell;
    const unsigned int n_verts_per_face = GeometryInfo<2>::vertices_per_cell;

    unsigned int n_vertices = mesh->get_n_vertices();
    unsigned int n_elems = mesh->get_n_hexahedra();
    unsigned int n_faces = mesh->get_n_quadrangles();
    unsigned int i;

    unsigned int material_id = 0;

    vector<Point<DIM> > vertices(n_vertices); // array for vertices
    vector<CellData<DIM> > cells(n_elems);    // array for elements
    SubCellData subcelldata;
    subcelldata.boundary_quads.reserve(n_faces); // array for faces

    for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
        for (i = 0; i < DIM; ++i)
            vertices[vertex](i) = mesh->get_vertex(vertex).get_coord(i);

    // copy faces
    for (unsigned int face = 0; face < n_faces; ++face) {
        subcelldata.boundary_quads[face] = CellData<2>();
        for (i = 0; i < n_verts_per_face; ++i)
            subcelldata.boundary_quads[face].vertices[i] = mesh->get_quadrangle(face).get_vertex(i);

//        // to make sure that the cast wont fail
//        Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
//                ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
//        // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
//        Assert(material_id < numbers::internal_face_boundary_id,
//                ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));
//
//        subcelldata.boundary_quads[face].boundary_id = static_cast<types::boundary_id>(material_id);
    }

    // copy hexahedra
    for (unsigned int elem = 0; elem < n_elems; ++elem) {
        cells[elem] = CellData<DIM>();
        for (i = 0; i < n_verts_per_elem; ++i)
            cells[elem].vertices[i] = mesh->get_hexahedron(elem).get_vertex(i);

//        // to make sure that the cast wont fail
//        Assert(material_id<= std::numeric_limits<types::material_id>::max(),
//                ExcIndexRange(material_id,0,std::numeric_limits<types::material_id>::max()));
//        // we use only material_ids in the range from 0 to numbers::invalid_material_id-1
//        Assert(material_id < numbers::invalid_material_id,
//                ExcIndexRange(material_id,0,numbers::invalid_material_id));
//
//        cells[elem].material_id = static_cast<types::material_id>(material_id);
    }

    // Check consistency of subcelldata
    Assert(subcelldata.check_consistency(DIM), ExcInternalError());

    // do some clean-up on vertices...
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    // ... and cells
    GridReordering<DIM, DIM>::invert_all_cells_of_negative_grid(vertices, cells);
//    GridReordering<DIM, DIM>::reorder_cells(cells);
    triangulation.create_triangulation_compatibility(vertices, cells, subcelldata);
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
    int n_elems = mesh->get_n_elems();
    triangulation.clear();

    if (n_elems < 1) return;

    for (i = 0; i < n_nodes_in_elem; ++i) {
        node = mesh->get_elem(0, i);
        for (j = 0; j < n_coords; ++j)
            vertices1[i](j) = mesh->get_node(node, j);
    }

    if (n_elems == 1) {
        GridGenerator::simplex(triangulation, vertices1);
        return;
    }

    for (i = 0; i < n_nodes_in_elem; ++i) {
        node = mesh->get_elem(1, i);
        for (j = 0; j < n_coords; ++j)
            vertices2[i](j) = mesh->get_node(node, j);
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
            node = mesh->get_elem(elem, i);
            for (j = 0; j < n_coords; ++j)
                vertices1[i](j) = mesh->get_node(node, j);
        }
        GridGenerator::simplex(tr1, vertices1);
        GridGenerator::merge_triangulations(tr1, triangulation, triangulation);
        tr1.clear();
    }
}

// Generate simple mesh for test purposes
void DealII::make_simple_mesh() {
    const unsigned int d = 3;
    Triangulation<DIM> tr1;
    Triangulation<DIM> tr2;

    vector<Point<DIM> > vertices1(DIM + 1);
    vector<Point<DIM> > vertices2(DIM + 1);
    vector<Point<DIM> > vertices3(DIM + 1);
    vector<Point<DIM> > vertices4(DIM + 1);

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
void DealII::mark_boundary(const AtomReader::Sizes* sizes, const AtomReader::Types* types) {
    typename Triangulation<DIM>::active_face_iterator face;

    // Loop through the faces and mark them according the location of its centre
    for (face = triangulation.begin_active_face(); face != triangulation.end(); ++face) {
        if (face->at_boundary()) {
            if (on_boundary(face->center()[0], sizes->xmin, sizes->xmax))
                face->set_all_boundary_ids(types->type_edge);
            else if (on_boundary(face->center()[1], sizes->ymin, sizes->ymax))
                face->set_all_boundary_ids(types->type_edge);
            else if (on_boundary(face->center()[2], sizes->zmaxbox, sizes->zmaxbox))
                face->set_all_boundary_ids(types->type_zmax);
            else
                face->set_all_boundary_ids(types->type_surf);
        }
    }
}

// Function to determine whether the center of face is on the boundary of simulation cell or not
bool DealII::on_boundary(const double face, const double face_min, const double face_max) {
    const double eps = 0.1;
    bool b1 = fabs(face - face_min) < eps;
    bool b2 = fabs(face_max - face) < eps;
    return b1 || b2;
}

// Insert boundary conditions to the system
void DealII::assemble_system(const AtomReader::Types* types) {
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
                cell_rhs(i) += fe_values.shape_value(i, q_index) * space_charge
                        * fe_values.JxW(q_index);

            // Assemble the matrix
            for (i = 0; i < dofs_per_cell; ++i)
                for (j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q_index)
                            * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index);

            // Cycle for faces of each cell
            for (unsigned int f = 0; f < faces_per_cell; ++f)
                // Check if face is located at top boundary  
                if (cell->face(f)->boundary_id() == types->type_zmax) {
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
    VectorTools::interpolate_boundary_values(dof_handler, types->type_surf, ZeroFunction<DIM>(),
            copper_boundary_value);

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
    gout.write_msh(triangulation, outfile);
}

} /* namespace femocs */
