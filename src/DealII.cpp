/*
 * DealII.cpp
 /*
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

// Modified read_msh function
void DealII::import_file_vol2(const string file_name) {
    ifstream in(file_name);
    AssertThrow(in, ExcIO());

    const unsigned int n_verts_per_elem = GeometryInfo<3>::vertices_per_cell;
    const unsigned int n_verts_per_face = GeometryInfo<2>::vertices_per_cell;

    unsigned int n_vertices;
    unsigned int n_cells;
    unsigned int n_elems;
    unsigned int n_faces;

    unsigned int dummy;
    unsigned int i;
    string line;

    in >> line;

    // first determine file format
    unsigned int gmsh_file_format;
    if (line == "$NOD")
        gmsh_file_format = 1;
    else if (line == "$MeshFormat")
        gmsh_file_format = 2;
    else
        gmsh_file_format = -1;

    Assert(gmsh_file_format == 2, ExcNotImplemented());

    double version;
    unsigned int file_type, data_size;

    in >> version >> file_type >> data_size;

    Assert((version >= 2.0) && (version <= 2.2), ExcNotImplemented());
    Assert(file_type == 0, ExcNotImplemented());
    Assert(data_size == sizeof(double), ExcNotImplemented());

    // read the end of the header and the first line of the nodes description to synch
    // ourselves with the format 1 handling above
    in >> line;

    in >> line;
    // if the next block is of kind $PhysicalNames, ignore it
    if (line == "$PhysicalNames") {
        do {
            in >> line;
        } while (line != "$EndPhysicalNames");
        in >> line;
    }

    // now read the nodes list
    in >> n_vertices;
    vector<Point<DIM> > vertices(n_vertices);
    // set up mapping between numbering in msh-file (nod) and in the vertices vector
    map<int, int> vertex_indices;
    for (unsigned int vertex = 0; vertex < n_vertices; ++vertex) {
        int vertex_number;
        double x[3];
        // read vertex
        in >> vertex_number >> x[0] >> x[1] >> x[2];
        for (i = 0; i < DIM; ++i)
            vertices[vertex](i) = x[i];
        // store mapping
        vertex_indices[vertex_number] = vertex;
    }

    // Assert we reached the end of the block
    in >> line;
    static const std::string end_nodes_marker[] = { "$ENDNOD", "$EndNodes" };
    AssertThrow(line == end_nodes_marker[gmsh_file_format - 1], ExcInternalError());

    // Now read in next bit
    in >> line;
    static const std::string begin_elements_marker[] = { "$ELM", "$Elements" };
    AssertThrow(line == begin_elements_marker[gmsh_file_format - 1], ExcInternalError());

    in >> n_cells;

    // set up array of cells and subcells (faces). In 1d, there is currently no
    // standard way in deal.II to pass boundary indicators attached to individual
    // vertices, so do this by hand via the boundary_ids_1d array
    vector<CellData<DIM> > cells;
    SubCellData subcelldata;

    for (unsigned int cell = 0; cell < n_cells; ++cell) {
        unsigned int cell_type;
        unsigned int material_id;

        in >> dummy >> cell_type;     // ELM-TYPE

        // read the tags; ignore all but the first one which we will
        // interpret as the material_id (for cells) or boundary_id
        // (for faces)
        unsigned int n_tags;
        in >> n_tags;
        if (n_tags > 0)
            in >> material_id;
        else
            material_id = 0;

        for (i = 1; i < n_tags; ++i)
            in >> dummy;

        // Hexahedron:
        if (cell_type == 5) {
            // allocate and read indices
            cells.push_back(CellData<DIM>());
            for (i = 0; i < n_verts_per_elem; ++i) {
                in >> cells.back().vertices[i];
            }


            // to make sure that the cast wont fail
            Assert(material_id<= std::numeric_limits<types::material_id>::max(),
                    ExcIndexRange(material_id,0,std::numeric_limits<types::material_id>::max()));
            // we use only material_ids in the range from 0 to numbers::invalid_material_id-1
            Assert(material_id < numbers::invalid_material_id,
                    ExcIndexRange(material_id,0,numbers::invalid_material_id));

            cells.back().material_id = static_cast<types::material_id>(material_id);

            // transform from ucd to consecutive numbering
            for (i = 0; i < n_verts_per_elem; ++i)
                cells.back().vertices[i] = vertex_indices[cells.back().vertices[i]];
        }

        // Quadrangle:
        else if (cell_type == 3) {
            // boundary info
            subcelldata.boundary_quads.push_back(CellData<2>());
            in >> subcelldata.boundary_quads.back().vertices[0]
                    >> subcelldata.boundary_quads.back().vertices[1]
                    >> subcelldata.boundary_quads.back().vertices[2]
                    >> subcelldata.boundary_quads.back().vertices[3];

            // to make sure that the cast wont fail
            Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
                    ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
            // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
            Assert(material_id < numbers::internal_face_boundary_id,
                    ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));

            subcelldata.boundary_quads.back().boundary_id =
                    static_cast<types::boundary_id>(material_id);

            // transform from gmsh to consecutive numbering
            for (i = 0; i < n_verts_per_face; ++i)
                if (vertex_indices.find(subcelldata.boundary_quads.back().vertices[i])
                        != vertex_indices.end())
                    // vertex with this index exists
                    subcelldata.boundary_quads.back().vertices[i] =
                            vertex_indices[subcelldata.boundary_quads.back().vertices[i]];
                else {
                    // no such vertex index
                    Assert(false,
                            ExcInvalidVertexIndex(cell, subcelldata.boundary_quads.back().vertices[i]));
                    subcelldata.boundary_quads.back().vertices[i] = numbers::invalid_unsigned_int;
                }
        }
    }

    // Assert we reached the end of the block
    in >> line;
    static const string end_elements_marker[] = { "$ENDELM", "$EndElements" };
    AssertThrow(line == end_elements_marker[gmsh_file_format - 1], ExcInternalError());

    // check that no forbidden arrays are used
    Assert(subcelldata.check_consistency(DIM), ExcInternalError());

    // do some clean-up on
    // vertices...
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    // ... and cells
//    GridReordering<DIM, DIM>::invert_all_cells_of_negative_grid(vertices, cells);
//    GridReordering<DIM, DIM>::reorder_cells(cells);
    triangulation.create_triangulation_compatibility(vertices, cells, subcelldata);
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
void DealII::mark_boundary() {
    typename Triangulation<DIM>::active_face_iterator face;

    // Loop through the faces and mark them according the location of its centre
    for (face = triangulation.begin_active_face(); face != triangulation.end(); ++face) {
        if (face->at_boundary()) {
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
                if (cell->face(f)->boundary_id() == types.top) {
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
    VectorTools::interpolate_boundary_values(dof_handler, types.surface, ZeroFunction<DIM>(),
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
