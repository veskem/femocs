

#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/work_stream.h>

#include "PoissonSolver.h"
#include "Globals.h"


namespace femocs {

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
        particles(NULL), conf(NULL), interpolator(NULL), applied_field(0), applied_potential(0)
        {}

template<int dim>
PoissonSolver<dim>::PoissonSolver(const ParticleSpecies* particles_,
        const Config::Field* conf_, const LinearHexahedra* interpolator_) :
        DealSolver<dim>(),
        particles(particles_), conf(conf_), interpolator(interpolator_),
        applied_field(0), applied_potential(0)
        {}

template<int dim>
void PoissonSolver<dim>::mark_mesh() {
    this->mark_boundary(BoundaryID::vacuum_top, BoundaryID::copper_surface,
            BoundaryID::vacuum_sides, BoundaryID::copper_surface);
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
void PoissonSolver<dim>::potential_efield_at(vector<double> &potentials, vector<Tensor<1, dim>> &fields,
        const vector<int> &cells, const vector<int> &verts) const
{
    this->solution_at(potentials, cells, verts);
    this->solution_grad_at(fields, cells, verts);
}

template<int dim>
void PoissonSolver<dim>::charge_dens_at(vector<double> &charge_dens,
        const vector<int> &cells, const vector<int> &verts)
{
    const int n_nodes = cells.size();
    require(n_nodes == verts.size(), "Invalid vectors sizes for cells and vertices: "
            + d2s(n_nodes) + ", " + d2s(verts.size()));

    charge_dens.resize(n_nodes);
    this->calc_dof_volumes();

    for (unsigned i = 0; i < n_nodes; i++) {
        // Using DoFAccessor (groups.google.com/forum/?hl=en-GB#!topic/dealii/azGWeZrIgR0)
        // NB: only works without refinement !!!
        typename DoFHandler<dim>::active_cell_iterator dof_cell(&this->triangulation,
                0, cells[i], &this->dof_handler);

        double rhs = this->system_rhs_save[dof_cell->vertex_dof_index(verts[i], 0)];
        double vol = this->dof_volume[dof_cell->vertex_dof_index(verts[i], 0)];
        charge_dens[i] = rhs / vol;
    }
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
void PoissonSolver<dim>::assemble(const bool full_run) {
    require(conf, "NULL conf can't be used!");
    require(conf->anode_BC == "neumann" || conf->anode_BC == "dirichlet",
            "Unimplemented anode BC: " + conf->anode_BC);

    if (full_run) this->system_matrix = 0;
    this->system_rhs = 0;

    if (conf->anode_BC == "neumann") {
        if (full_run) {
            assemble_parallel();
            this->append_dirichlet(BoundaryID::copper_surface, this->dirichlet_bc_value);
        }
        this->assemble_rhs(BoundaryID::vacuum_top);
    } else {
        if (full_run) {
            assemble_parallel();
            this->append_dirichlet(BoundaryID::copper_surface, this->dirichlet_bc_value);
            this->append_dirichlet(BoundaryID::vacuum_top, applied_potential);
        }
    }

    if (conf->mode != "laplace") assemble_space_charge();
    // save the system right-hand-side to make it possible to write space charge into file
    if (this->write_time()) this->system_rhs_save = this->system_rhs;
    if (full_run) this->apply_dirichlet();
}

template<int dim>
void PoissonSolver<dim>::assemble_parallel() {
    LinearSystem system(&this->system_rhs, &this->system_matrix);
    QGauss<dim> quadrature_formula(this->quadrature_degree);

    const unsigned int n_dofs = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    WorkStream::run(this->dof_handler.begin_active(),this->dof_handler.end(),
            std_cxx11::bind(&PoissonSolver<dim>::assemble_local_cell,
                    this,
                    std_cxx11::_1,
                    std_cxx11::_2,
                    std_cxx11::_3),
            std_cxx11::bind(&PoissonSolver<dim>::copy_global_cell,
                    this,
                    std_cxx11::_1,
                    std_cxx11::ref(system)),
            ScratchData(this->fe, quadrature_formula, update_gradients | update_quadrature_points | update_JxW_values),
            CopyData(n_dofs, n_q_points)
    );
}

template<int dim>
void PoissonSolver<dim>::assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
        ScratchData &scratch_data, CopyData &copy_data) const
{
    const unsigned int n_dofs = copy_data.n_dofs;
    const unsigned int n_q_points = copy_data.n_q_points;

    scratch_data.fe_values.reinit(cell);

    // Local matrix assembly
    copy_data.cell_matrix = 0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
        for (unsigned int i = 0; i < n_dofs; ++i) {
            for (unsigned int j = 0; j < n_dofs; ++j) {
                copy_data.cell_matrix(i, j) += scratch_data.fe_values.JxW(q) *
                scratch_data.fe_values.shape_grad(i, q) * scratch_data.fe_values.shape_grad(j, q);
            }
        }
    }

    // Nothing to add to local right-hand-side vector assembly
//    copy_data.cell_rhs = 0;

    // Obtain dof indices for updating global matrix and right-hand-side vector
    cell->get_dof_indices(copy_data.dof_indices);
}

// In general case, use only Deal.II built-in tools
template<int dim>
void PoissonSolver<dim>::assemble_space_charge() {
    if (!particles) {
        write_silent_msg("No charged particles present!");
        return;
    }

    if (dim == 3 && this->shape_degree == 1) {
        assemble_space_charge_fast();
        return;
    }

    vector<types::global_dof_index> local_dof_indices(this->fe.dofs_per_cell);

    // loop over particles
    for (SuperParticle const &particle : *particles) {
        Point<dim> p_deal = Point<dim>(particle.pos.x, particle.pos.y, particle.pos.z);
        //get particle's active cell iterator
        typename DoFHandler<dim>::active_cell_iterator cell(&this->triangulation, 0, particle.cell, &this->dof_handler);

        //get the node indices of the particle's cell
        cell->get_dof_indices(local_dof_indices);

        //get the shape functions of the cell on the given point
        vector<double> sf = this->shape_funs(p_deal, particle.cell);

        //loop over nodes of the cell and add the particle's charge to the system rhs
        for (unsigned int i = 0; i < this->fe.dofs_per_cell; ++i)
            this->system_rhs(local_dof_indices[i]) += sf[i] * particles->q_over_eps0 * particles->get_Wsp();
    }
}

// In linear 3D case, Femocs interpolator can be used to speed up calculations
template<>
void PoissonSolver<3>::assemble_space_charge_fast() {
    static constexpr int n_dofs = GeometryInfo<3>::vertices_per_cell;
    const double charge_factor = particles->q_over_eps0 * particles->get_Wsp();
    vector<types::global_dof_index> local_dof_indices(n_dofs);

    // loop over particles
    for (SuperParticle const &particle : *particles) {
        // get the node indices of the particle's cell
        typename DoFHandler<3>::active_cell_iterator cell(&this->triangulation, 0, particle.cell, &this->dof_handler);
        cell->get_dof_indices(local_dof_indices);

        //get the shape functions of the cell on the given point
        int femocs_cell = interpolator->deal2femocs(particle.cell);
        array<double,n_dofs> shape_fun = interpolator->shape_funs_dealii(particle.pos, femocs_cell);

        //loop over nodes of the cell and add the particle's charge to the system rhs
        for (int i = 0; i < n_dofs; ++i)
            this->system_rhs(local_dof_indices[i]) += shape_fun[i] * charge_factor;
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

} // namespace femocs
