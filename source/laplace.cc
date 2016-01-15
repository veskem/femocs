
#include "laplace.h"
#include "utility.h"


namespace Emitter {

Laplace::Laplace() :
		fe(2), dof_handler(triangulation) {
}


class LaplacePostProcessor : public DataPostprocessorVector<2> {
public:
	LaplacePostProcessor();

	virtual void
    compute_derived_quantities_scalar (const std::vector<double>            	&uh,
                                       const std::vector<Tensor<1,2> >  		&duh,
                                       const std::vector<Tensor<2,2> >  		&dduh,
                                       const std::vector<Point<2> >        		&normals,
                                       const std::vector<Point<2> >          	&evaluation_points,
                                       std::vector<Vector<double> >          	&computed_quantities) const;
};

LaplacePostProcessor::LaplacePostProcessor() : DataPostprocessorVector<2>("Field", update_values | update_gradients) {}

void LaplacePostProcessor::
compute_derived_quantities_scalar (	const std::vector<double>            	&uh,
									const std::vector<Tensor<1,2> >  		&duh,
									const std::vector<Tensor<2,2> >  		&/*dduh*/,
									const std::vector<Point<2> >        	&/*normals*/,
									const std::vector<Point<2> >          	&/*evaluation_points*/,
									std::vector<Vector<double> >          	&computed_quantities) const {
	Assert(computed_quantities.size() == uh.size(),
	       ExcDimensionMismatch (computed_quantities.size(), uh.size()));
	for (unsigned int i=0; i<computed_quantities.size(); i++) {
	    Assert(computed_quantities[i].size() == 2,
	           ExcDimensionMismatch (computed_quantities[i].size(), 2));
	    for (unsigned int d=0; d<2; ++d)
	    	computed_quantities[i](d) = duh[i][d];

	}
}

// Makes a geometry with protrusion-like shape
// NB: This is old version and boundary indexes might be incorrect
void Laplace::make_grid_old() {
	// Create a 3x3 grid
	GridGenerator::subdivided_hyper_cube(triangulation, 3, -1, 1);

	// Cell removal (only the middle bottom one)
	std::set<Triangulation<2>::active_cell_iterator> removal;
	Triangulation<2>::active_cell_iterator cell;
	for (cell = triangulation.begin_active(); cell != triangulation.end();
			++cell) {
		Point<2> v = cell->center();

		// If it is the bottom center cell, remove it
		if (v(0) > -0.5 && v(0) < 0.5 && v(1) < -0.1) {
			//std::cout << v(0) << " " << v(1) << std::endl;
			removal.insert(cell);
		}
	}
	// Remove the marked cells
	GridGenerator::create_triangulation_with_removed_cells(triangulation,
			removal, triangulation);

	// Loop over cells, define manifolds and boundaries
	for (cell = triangulation.begin_active(); cell != triangulation.end();
			++cell) {
		Point<2> v = cell->center();

		// If it is the middle cell, define the spherical manifold on it
		if (v(0) > -0.5 && v(0) < 0.5 && v(1) > -0.1 && v(1) < 0.1) {
			cell->set_all_manifold_ids(10);
		}

		if (v(0) < -0.1)
			cell->face(0)->set_boundary_id(1);
		if (v(0) > 0.1)
			cell->face(1)->set_boundary_id(1);
		if (v(1) > 0.1)
			cell->face(3)->set_boundary_id(2);
	}
	// Define the spherical manifold
	const SphericalManifold<2> manifold_desc(Point<2>(0, -0.334));
	triangulation.set_manifold(10, manifold_desc);

	triangulation.refine_global(3);
	triangulation.execute_coarsening_and_refinement();

	output_mesh(triangulation, "mesh.eps");
	triangulation.set_manifold(10);
}

void Laplace::make_grid() {

	double radius = 0.99;
	double outer_radius = 9.0;
	double h = 4.0;

	double eps = 1e-6;

	unsigned char tip_manifold = 10;


	Triangulation<2> mesh;

	// --------------------------------------------------------
	// Mesh of the curved tip
	// --------------------------------------------------------

	// The inner circle has boundary id of 1 and outer square has 0 (colorize=false)
	GridGenerator::hyper_cube_with_cylindrical_hole(mesh, radius, outer_radius,
			0.0, 0, false);

	// Remove half of the cells
	std::set<Triangulation<2>::active_cell_iterator> removal;
	Triangulation<2>::active_cell_iterator cell;
	for (cell = mesh.begin_active(); cell != mesh.end(); ++cell) {
		Point<2> v = cell->center();
		// remove all cells, which are under line y=0
		if (v(1) < 0) {
			removal.insert(cell);
		}
	}

	// Remove the marked cells
	// NB: This removes the initial boundary_id's generated with GridGenerator
	GridGenerator::create_triangulation_with_removed_cells(mesh, removal, mesh);

	// --------------------------------------------------------
	// Mesh surrounding the curved tip
	// --------------------------------------------------------

	// bottom left and right rectangles
	Triangulation<2> mesh_h;
	GridGenerator::hyper_rectangle(mesh_h, Point<2>(-outer_radius, -h),
			Point<2>(-radius, 0.0));
	GridGenerator::merge_triangulations(mesh, mesh_h, mesh);

	GridTools::shift(Point<2>(outer_radius + radius, 0.0), mesh_h);
	GridGenerator::merge_triangulations(mesh, mesh_h, triangulation);


	// --------------------------------------------------------
	// Define boundary IDs and manifolds
	// --------------------------------------------------------

	for (cell = triangulation.begin_active(); cell != triangulation.end(); ++cell) {
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary()) {
				Point<2> center = cell->face(f)->center();

				// The face is on the inner circle
				if (center.norm() - eps < radius) {
					cell->face(f)->set_all_boundary_ids(copper_boundary);
					cell->face(f)->set_all_manifold_ids(tip_manifold);
				}

				// The face is on bottom or on the side of the protrusion
				if (center[1] - eps < -h
						|| (center[1] - eps < 0
								&& std::abs(center[0]) - eps < radius))
					cell->face(f)->set_all_boundary_ids(copper_boundary);

				// The face is on top
				if (center[1] + eps > outer_radius)
					cell->face(f)->set_all_boundary_ids(vacuum_top);
			}
		}
	}

	const SphericalManifold<2> manifold_desc(Point<2>(0.0, 0.0));
	triangulation.set_manifold(tip_manifold, manifold_desc);


	triangulation.refine_global(5);
	//mesh.execute_coarsening_and_refinement();

	output_mesh(triangulation, "laplace_mesh.eps");

	// before using the mesh further, all manifolds must be reset to FlatManifold
	triangulation.set_manifold(tip_manifold);

}

void Laplace::setup_system() {
	dof_handler.distribute_dofs(fe);
	std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
			<< std::endl;

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);

	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

void Laplace::assemble_system() {
	QGauss<2> quadrature_formula(2);
	QGauss<1> face_quadrature_formula(2);

	FEValues<2> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values);

	FEFaceValues<2> fe_face_values(fe, face_quadrature_formula,
			update_values | update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) {
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
					cell_matrix(i, j) += (fe_values.shape_grad(i, q_index)
							* fe_values.shape_grad(j, q_index)
							* fe_values.JxW(q_index));

				// No space charge, so f=0
				cell_rhs(i) += (fe_values.shape_value(i, q_index) * 0
						* fe_values.JxW(q_index));
			}
		}

		// Neumann boundary condition 1.0 V/nm at faces with boundary_id vacuum_top
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary()
					&& cell->face(f)->boundary_id() == vacuum_top) {
				fe_face_values.reinit(cell, f);

				for (unsigned int q_index = 0; q_index < n_face_q_points;
						++q_index) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						cell_rhs(i) += (fe_face_values.shape_value(i, q_index)
								* 1.0 * fe_face_values.JxW(q_index));
					}
				}
			}

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	// Setting boundary values //

	// 0 potential at the "copper surface"
	std::map<types::global_dof_index, double> copper_boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, copper_boundary, ZeroFunction<2>(),
			copper_boundary_values);

	MatrixTools::apply_boundary_values(copper_boundary_values, system_matrix,
			solution, system_rhs);

	// The following code adds a Dirichlet boundary condition to boundaries with id 2,
	// but as Neumann conditions are used there, this is ignored.

//  // 1 V at the vacuum boundary
//  std::map<types::global_dof_index,double> vacuum_boundary_values;
//  VectorTools::interpolate_boundary_values (dof_handler,
//                                            vacuum_top,
//                                            ConstantFunction<2>(1),
//                                            vacuum_boundary_values);
//  MatrixTools::apply_boundary_values (vacuum_boundary_values,
//                                      system_matrix,
//                                      solution,
//                                      system_rhs);
}

void Laplace::solve() {
	SolverControl solver_control(3000, 1e-12);
	SolverCG<> solver(solver_control);

	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

double Laplace::probe_field(const Point<2> &p) const {
	return VectorTools::point_gradient (dof_handler, solution, p).norm();
}

void Laplace::output_results() const {
	LaplacePostProcessor field_calculator; // needs to be before data_out
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "potential_v");
	data_out.add_data_vector(solution, field_calculator);
	data_out.build_patches();

	Point<2> eval_point(0.0, 1.0);

	std::cout << "Sol  at (" << eval_point << "): "
			<< VectorTools::point_value (dof_handler, solution, eval_point)
			<< std::endl;
	std::cout << "Grad at (" << eval_point << "): "
			<< probe_field(eval_point)
			<< std::endl;

	std::ofstream output("field_solution.vtk");
	data_out.write_vtk(output);
}

void Laplace::run() {
	make_grid();
	setup_system();
	assemble_system();
	solve();
	output_results();
}

} //end namespace
