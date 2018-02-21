

#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>

#include "PoissonSolver.h"
#include "Globals.h"



namespace fch {
using namespace dealii;
using namespace std;

// ----------------------------------------------------------------------------------------
/* Class for outputting the field distribution
 * being calculated from potential distribution */
template <int dim>
class LaplacePostProcessor : public DataPostprocessorVector<dim> {
public:
    LaplacePostProcessor() : DataPostprocessorVector<dim>("field", update_gradients) {}

    void
    compute_derived_quantities_scalar (
            const vector<double>             &/*uh*/,
            const vector<Tensor<1,dim> >     &duh,
            const vector<Tensor<2,dim> >     &/*dduh*/,
            const vector<Point<dim> >        &/*normals*/,
            const vector<Point<dim> >        &/*evaluation_points*/,
            vector<Vector<double> >          &computed_quantities) const {
        for (unsigned int i=0; i<computed_quantities.size(); i++) {
            for (unsigned int d=0; d<dim; ++d)
                computed_quantities[i](d) = duh[i][d];

        }
    }
};
// ----------------------------------------------------------------------------

template<int dim>
PoissonSolver<dim>::PoissonSolver() : DealSolver<dim>(),
        particles(NULL), conf(NULL), applied_field(0), applied_potential(0)
        {};

template<int dim>
PoissonSolver<dim>::PoissonSolver(const ParticleSpecies* particles_, const femocs::Config::Field *conf_) :
        DealSolver<dim>(),
        particles(particles_), conf(conf_), applied_field(0), applied_potential(0)
        {};

template<int dim>
void PoissonSolver<dim>::mark_mesh() {
//    MeshPreparer<dim> mesh_preparer;
//    mesh_preparer.mark_vacuum_boundary(&this->triangulation);
    this->mark_boundary(BoundaryId::vacuum_top, BoundaryId::copper_surface,
            BoundaryId::vacuum_sides, BoundaryId::copper_surface);
}

template<int dim>
double PoissonSolver<dim>::probe_efield_norm(const Point<dim> &p) const {
    return VectorTools::point_gradient (this->dof_handler, this->solution, p).norm();
}

template<int dim>
double PoissonSolver<dim>::probe_potential(const Point<dim> &p) const{
    return VectorTools::point_value(this->dof_handler, this->solution, p);
}

template<int dim>
double PoissonSolver<dim>::probe_potential(const Point<dim> &p, int cell_index) const {
    return probe_potential(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
double PoissonSolver<dim>::probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const {

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&this->triangulation, 0, max(0,cell_index), &this->dof_handler);

    //point in transformed unit cell coordinates
    Point<dim> p_cell;

    if (cell_index < 0){ // in case the cell index is unknown (argument cell_index < 0)
        const pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, this->dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    }else // cell index is known
        p_cell = mapping.transform_real_to_unit_cell(cell, p);


    //create virtual quadrature point
    const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));

    FEValues<dim> fe_values(mapping, this->fe, quadrature, update_values);
    fe_values.reinit(cell);

    vector<Vector<double> > u_value(1, Vector<double> (this->fe.n_components()));
    fe_values.get_function_values(this->solution, u_value);

    return u_value[0][0];
}

template<int dim>
double PoissonSolver<dim>::probe_efield_norm(const Point<dim> &p, int cell_index) const {
    return probe_efield(p, cell_index, StaticMappingQ1<dim,dim>::mapping).norm();
}

template<int dim>
Tensor<1, dim, double> PoissonSolver<dim>::probe_efield(const Point<dim> &p, int cell_index) const {
    return probe_efield(p, cell_index, StaticMappingQ1<dim,dim>::mapping);
}

template<int dim>
Tensor<1, dim, double> PoissonSolver<dim>::probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const {

    static double tconstruct = 0, tcalc = 0, t0;

    //get active cell iterator from cell index
    typename DoFHandler<dim>::active_cell_iterator cell(&this->triangulation, 0, max(0,cell_index), &this->dof_handler);

    // transform the point from real to unit cell coordinates
    Point<dim> p_cell;
    if (cell_index < 0){
        const pair<typename DoFHandler<dim,dim>::active_cell_iterator, Point<dim> > cell_point
        = GridTools::find_active_cell_around_point (mapping, this->dof_handler, p);
        cell = cell_point.first;
        p_cell = cell_point.second;
    }else
        p_cell = mapping.transform_real_to_unit_cell(cell, p);

    const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));

    FEValues<dim> fe_values(mapping, this->fe, quadrature, update_gradients);
    fe_values.reinit(cell);

    vector<vector<Tensor<1, dim, double> > >
    u_gradient(1, vector<Tensor<1, dim, double> > (this->fe.n_components()));
    fe_values.get_function_gradients(this->solution, u_gradient);

    return -u_gradient[0][0];
}

template<int dim>
vector<double> PoissonSolver<dim>::get_potential(const vector<int> &cell_indexes,
        const vector<int> &vert_indexes) {

    // Initialise potentials with a value that is immediately visible if it's not changed to proper one
    vector<double> potentials(cell_indexes.size(), 1e15);

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&this->triangulation,
                0, cell_indexes[i], &this->dof_handler);

        double potential = this->solution[dof_cell->vertex_dof_index(vert_indexes[i], 0)];
        potentials[i] = potential;
    }
    return potentials;
}

template<int dim>
vector<Tensor<1, dim> > PoissonSolver<dim>::get_efield(const vector<int> &cell_indexes,
        const vector<int> &vert_indexes) const {

    QGauss<dim> quadrature_formula(this->quadrature_degree);
    FEValues<dim> fe_values(this->fe, quadrature_formula, update_gradients);

    vector< Tensor<1, dim> > solution_gradients (quadrature_formula.size());

    vector<Tensor<1, dim> > fields(cell_indexes.size());

    for (unsigned i = 0; i < cell_indexes.size(); i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&this->triangulation,
                0, cell_indexes[i], &this->dof_handler);

        fe_values.reinit(dof_cell);
        fe_values.get_function_gradients(this->solution, solution_gradients);
        Tensor<1, dim> field = -1.0 * solution_gradients.at(vert_indexes[i]);

        fields[i] = field;
    }
    return fields;
}

template<int dim>
double PoissonSolver<dim>::get_face_bc(const unsigned int face) const {
    return applied_field;
}

template<int dim>
void PoissonSolver<dim>::setup(const double field, const double potential) {
    DealSolver<dim>::setup_system();
    applied_field = field;
    applied_potential = potential;
}

template<int dim>
void PoissonSolver<dim>::assemble_laplace(const bool first_time) {
    require(conf, "NULL conf can't be used!");

    this->system_rhs = 0;

    if (first_time) {
        assemble_lhs();
        this->append_dirichlet(fch::BoundaryId::copper_surface, 0.);
    } else
        this->restore_system();

    this->assemble_rhs(fch::BoundaryId::vacuum_top);
    this->apply_dirichlet();
}

template<int dim>
void PoissonSolver<dim>::assemble_poisson(const bool first_time) {
    require(conf, "NULL conf can't be used!");
    require(conf->anodeBC == "neumann" || conf->anodeBC == "dirichlet",
            "Unimplemented anode BC: " + conf->anodeBC);

    this->system_rhs = 0;

    if (conf->anodeBC == "neumann") {
        if (first_time) {
            assemble_lhs();
            this->append_dirichlet(fch::BoundaryId::copper_surface, this->dirichlet_bc_value);
        } else
            this->restore_system();
        this->assemble_rhs(fch::BoundaryId::vacuum_top);

    } else {
        if (first_time) {
            assemble_lhs();
            this->append_dirichlet(fch::BoundaryId::copper_surface, this->dirichlet_bc_value);
            this->append_dirichlet(fch::BoundaryId::vacuum_top, applied_potential);
        } else
            this->restore_system();
    }

    assemble_space_charge();
    this->apply_dirichlet();
}

template<int dim>
void PoissonSolver<dim>::assemble_lhs() {
    this->system_matrix = 0;

    QGauss<dim> quadrature_formula(this->quadrature_degree);
    FEValues<dim> fe_values(this->fe, quadrature_formula,
            update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Iterate over all cells (quadrangles in 2D, hexahedra in 3D) of the mesh
    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active();
    for (; cell != this->dof_handler.end(); ++cell) {
        fe_values.reinit(cell);

        // Assemble system matrix elements corresponding the current cell
        cell_matrix = 0;
        for (unsigned int q = 0; q < n_q_points; ++q) {
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
                            * fe_values.JxW(q);
            }
        }

        // Add the current cell matrix and rhs entries to the system sparse matrix
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                this->system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
    
    this->save_system();
}

template<int dim>
void PoissonSolver<dim>::assemble_space_charge() {
    if (!particles) {
        femocs::write_verbose_msg("No charged particles present!");
        return;
    }

    vector<types::global_dof_index> local_dof_indices(this->fe.dofs_per_cell);

    // loop over particles
    for (auto particle : particles->parts) {
        Point<dim> p_deal = Point<dim>(particle.pos.x, particle.pos.y, particle.pos.z);
        //get particle's active cell iterator
        typename DoFHandler<dim>::active_cell_iterator cell(&this->triangulation, 0, particle.cell, &this->dof_handler);

        //get the node indices of the particle's cell
        cell->get_dof_indices(local_dof_indices);

        //get the shape functions of the cell on the given point
        vector<double> sf = this->shape_funs(p_deal, particle.cell);

        //loop over nodes of the cell and add the particle's charge to the system rhs
        for (int i = 0; i < this->fe.dofs_per_cell; ++i)
            this->system_rhs(local_dof_indices[i]) += sf[i] * particles->q_over_eps0 * particles->Wsp;
    }
}

template<int dim>
void PoissonSolver<dim>::write_vtk(ofstream& out) const {
    LaplacePostProcessor<dim> post_processor; // needs to be before data_out
    DataOut<dim> data_out;

    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "electric_potential");
    data_out.add_data_vector(this->solution, post_processor);

    data_out.build_patches();
    data_out.write_vtk(out);
}

template class PoissonSolver<3> ;

} // namespace fch
