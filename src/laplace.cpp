


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

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
#include <deal.II/lac/precondition.h>
#include <deal.II/fe/fe_update_flags.h>


#include "laplace.h"


namespace fch {
using namespace dealii;




// ----------------------------------------------------------------------------------------
// Class for outputting the resulting field distribution (calculated from potential distr.)

template <int dim>
class LaplacePostProcessor : public DataPostprocessorVector<dim> {
public:
    LaplacePostProcessor() : DataPostprocessorVector<dim>("Field", update_values | update_gradients) {}

    void
    compute_derived_quantities_scalar ( const std::vector<double>             &/*uh*/,
            const std::vector<Tensor<1,dim> >     &duh,
            const std::vector<Tensor<2,dim> >     &/*dduh*/,
            const std::vector<Point<dim> >          &/*normals*/,
            const std::vector<Point<dim> >          &/*evaluation_points*/,
            std::vector<Vector<double> >            &computed_quantities) const {
        for (unsigned int i=0; i<computed_quantities.size(); i++) {
            for (unsigned int d=0; d<dim; ++d)
                computed_quantities[i](d) = duh[i][d];

        }
    }
};
// ----------------------------------------------------------------------------


template<int dim>
Laplace<dim>::Laplace() : fe(shape_degree), dof_handler(triangulation) {}

template<int dim>
Triangulation<dim>* Laplace<dim>::get_triangulation() {
    return &triangulation;
}

template <int dim>
DoFHandler<dim>* Laplace<dim>::get_dof_handler() {
    return &dof_handler;
}


template<int dim>
void Laplace<dim>::import_mesh_from_file(const std::string file_name) {
    MeshPreparer<dim> mesh_preparer;

    mesh_preparer.import_mesh_from_file(&triangulation, file_name);
    mesh_preparer.mark_vacuum_boundary(&triangulation);
}

template<int dim>
bool Laplace<dim>::import_mesh_directly(std::vector<Point<dim> > vertices, std::vector<CellData<dim> > cells) {

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
    mesh_preparer.mark_vacuum_boundary(&triangulation);

    return true;
}

template<int dim>
void Laplace<dim>::output_mesh(const std::string file_name) {
    MeshPreparer<dim> mesh_preparer;
    mesh_preparer.output_mesh(&triangulation, file_name);
}

template<int dim>
double Laplace<dim>::probe_efield_norm(const Point<dim> &p) const {
    return VectorTools::point_gradient (dof_handler, solution, p).norm();
}

template<int dim>
double Laplace<dim>::probe_potential(const Point<dim> &p) const{
    return VectorTools::point_value(dof_handler, solution, p);
}

template<int dim>
double Laplace<dim>::probe_potential(const Point<dim> &p, int cell_index) const {
    return probe_potential(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
double Laplace<dim>::probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const {

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, max(0,cell_index), &dof_handler);

    //point in transformed unit cell coordinates
    Point<dim> p_cell;

    if (cell_index < 0){ // in case the cell index is unknown (argument cell_index < 0)
        const std::pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    }else // cell index is known
        p_cell = mapping.transform_real_to_unit_cell(cell, p);


    //create virtual quadrature point
    const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));

    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell);

    std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
    fe_values.get_function_values(solution, u_value);

    return u_value[0][0];
}


template<int dim>
double Laplace<dim>::probe_efield_norm(const Point<dim> &p, int cell_index) const {
    return probe_efield(p, cell_index, StaticMappingQ1<dim,dim>::mapping).norm();
}

template<int dim>
Tensor<1, dim, double> Laplace<dim>::probe_efield(const Point<dim> &p, int cell_index) const {
    return probe_efield(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
Tensor<1, dim, double> Laplace<dim>::probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const {

    static double tconstruct = 0, tcalc = 0, t0;

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, max(0,cell_index), &dof_handler);

    // transform the point from real to unit cell coordinates
    Point<dim> p_cell;
    if (cell_index < 0){
        const std::pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    }else
        p_cell = mapping.transform_real_to_unit_cell(cell, p);

    const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));

    FEValues<dim> fe_values(mapping, fe, quadrature, update_gradients);
    fe_values.reinit(cell);

    std::vector<std::vector<Tensor<1, dim, double> > >
    u_gradient(1, std::vector<Tensor<1, dim, double> > (fe.n_components()));
    fe_values.get_function_gradients(solution, u_gradient);

    return -u_gradient[0][0];
}

template<int dim>
std::vector<double> Laplace<dim>::shape_funs(const Point<dim> &p, int cell_index) const {
    return shape_funs(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}


template<int dim>
std::vector<double> Laplace<dim>::shape_funs(const Point<dim> &p, const int cell_index,
                                               Mapping<dim,dim>& mapping) const {

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, max(0,cell_index), &dof_handler);
    //point in transformed unit cell coordinates
    Point<dim> p_cell;

    if (cell_index < 0){ // in case the cell index is unknown (argument cell_index < 0)
        const std::pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    }else // cell index is known
        p_cell = mapping.transform_real_to_unit_cell(cell, p);

    Point<dim> p_unit_cell = GeometryInfo<dim>::project_to_unit_cell(p_cell);

    //create virtual quadrature point
    const Quadrature<dim> quadrature(p_unit_cell);

    //define fevalues object
    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);

    fe_values.reinit(cell);

    std::vector<double> sfuns(fe.dofs_per_cell);

    for (int i = 0; i < sfuns.size(); i++){
        sfuns[i] = fe_values.shape_value(i,0);
    }

    return sfuns;
}

template<int dim>
double Laplace<dim>::get_cell_vol(int cellid){

    typename DoFHandler<dim>::active_cell_iterator cell(&triangulation,
            0, cellid, &dof_handler);
    return cell->measure();
}

template<int dim>
std::vector<double> Laplace<dim>::get_potential(const std::vector<int> &cell_indexes,
        const std::vector<int> &vert_indexes) {

    // Initialise potentials with a value that is immediately visible if it's not changed to proper one
    std::vector<double> potentials(cell_indexes.size(), 1e15);

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation,
                0, cell_indexes[i], &dof_handler);

        double potential = solution[dof_cell->vertex_dof_index(vert_indexes[i], 0)];
        potentials[i] = potential;
    }
    return potentials;
}

template<int dim>
std::vector<Tensor<1, dim> > Laplace<dim>::get_efield(const std::vector<int> &cell_indexes,
        const std::vector<int> &vert_indexes) const {

    QGauss<dim> quadrature_formula(quadrature_degree);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients);

    std::vector< Tensor<1, dim> > solution_gradients (quadrature_formula.size());

    std::vector<Tensor<1, dim> > fields(cell_indexes.size());

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&triangulation,
                0, cell_indexes[i], &dof_handler);

        fe_values.reinit(dof_cell);
        fe_values.get_function_gradients(solution, solution_gradients);
        Tensor<1, dim> field = -1.0 * solution_gradients.at(vert_indexes[i]);

        fields[i] = field;
    }
    return fields;
}

template<int dim>
void Laplace<dim>::setup_system(bool first_time) {

    if (first_time){ // find n_dofs
        dof_handler.distribute_dofs(fe);
    }

    system_rhs.reinit(dof_handler.n_dofs()); // set rhs to zeros

    if (!first_time){
        return;
    }

    boundary_values.clear();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    system_matrix_save.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
}

template<int dim>
void Laplace<dim>::assemble_system_lhs() {

    QGauss<dim> quadrature_formula(quadrature_degree);

    FEValues<dim> fe_values(fe, quadrature_formula,
            update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    // Iterate over all cells (quadrangles in 2D, hexahedra in 3D) of the mesh
    for (; cell != endc; ++cell) {
        fe_values.reinit(cell);
        cell_matrix = 0;

        // Assemble system matrix elements corresponding the current cell
        for (unsigned int q = 0; q < n_q_points; ++q) {
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
                            * fe_values.JxW(q));

            }
        }

        // Add the current cell matrix and rhs entries to the system sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
}

template<int dim>
void Laplace<dim>::assemble_system_neuman(BoundaryId bid, double applied_field) {

    QGauss<dim-1> face_quadrature_formula(quadrature_degree);

    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
            update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    // Iterate over all cells (quadrangles in 2D, hexahedra in 3D) of the mesh
    for (; cell != endc; ++cell) {
        cell_rhs = 0;

        // Apply Neumann boundary condition at faces on top of vacuum domain
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == bid) {
                fe_face_values.reinit(cell, f);

                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        cell_rhs(i) += (fe_face_values.shape_value(i, q)
                                * applied_field * fe_face_values.JxW(q));
                    }
                }
            }
        }

        // Add the current cell matrix and rhs entries to the system sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
}

template<int dim>
void Laplace<dim>::assemble_system_pointcharge(femocs::ParticleSpecies &particles) {

    std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
//    std::cout << "particles.Wsp = " << particles.get_Wsp() << std::endl;

    for(auto particle : particles.parts){ // loop over particles
        Point<dim> p_deal = Point<dim>(particle.pos.x, particle.pos.y, particle.pos.z);
        //get particle's active cell iterator
        typename DoFHandler<dim>::active_cell_iterator cell(&triangulation, 0, particle.cell, &dof_handler);

        //get the node indices of the particle's cell
        cell->get_dof_indices(local_dof_indices);

        //get the shape functions of the cell on the given point
        std::vector<double> sf = shape_funs(p_deal, particle.cell);

        //loop over nodes of the cell and add the particle's charge to the system rhs
        for (int j = 0; j < fe.dofs_per_cell; ++j){
            system_rhs(local_dof_indices[j]) += sf[j] * particles.q_over_eps0 * particles.get_Wsp();
        }

    }
}

template<int dim>
void Laplace<dim>::assemble_system_dirichlet(BoundaryId bid, double potential) {
    VectorTools::interpolate_boundary_values(dof_handler, bid, ConstantFunction<dim>(potential), boundary_values);

}

template<int dim>
int Laplace<dim>::solve(int max_iter, double tol, bool pc_ssor, double ssor_param) {

    SolverControl solver_control(max_iter, tol);
    SolverCG<> solver(solver_control);

    if (pc_ssor) {
        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(system_matrix, ssor_param);
        solver.solve(system_matrix, solution, system_rhs, preconditioner);
    } else {
        solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    }
    return solver_control.last_step();
}

template<int dim>
void Laplace<dim>::output_results(const std::string filename) const {
    LaplacePostProcessor<dim> field_calculator; // needs to be before data_out
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "potential_v");
    data_out.add_data_vector(solution, field_calculator);

    data_out.build_patches();

    try {
        std::ofstream output(filename);
        data_out.write_vtk(output);
    } catch (...) {
        std::cerr << "WARNING: Couldn't open " + filename << ". ";
        std::cerr << "Output is not saved." << std::endl;
    }
}

//template class Laplace<2> ;
template class Laplace<3> ;

} // namespace fch
