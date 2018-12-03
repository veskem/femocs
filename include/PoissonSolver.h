/*
 * laplace.h -> Laplace.h
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan, Mihkel
 */

#ifndef LAPLACE_H_
#define LAPLACE_H_

#include "DealSolver.h"
#include "Config.h"
#include "InterpolatorCells.h"
#include "ParticleSpecies.h"

namespace femocs {
using namespace dealii;
using namespace std;

/** @brief Class to solve Laplace equation in 2D or 3D
 * It is inspired by the step-3 of Deal.II tutorial
 * https://www.dealii.org/8.5.0/doxygen/deal.II/step_3.html
 */
template<int dim>
class PoissonSolver : public DealSolver<dim> {
public:
    PoissonSolver();
    PoissonSolver(const ParticleSpecies* particles, const Config::Field* conf, const LinearHexahedra* interpolator);

    /** Change the pointer to charged particles */
    void set_particles(const ParticleSpecies* parts) { particles = parts; }

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p) const;

    /** get the electric field norm at the specified point using dealii
     * (slow as it looks for the surrounding cell) */
    double probe_efield_norm(const Point<dim> &p, int cell_index) const;

    /** get the potential value at a specified point using dealii (slow)  */
    double probe_potential(const Point<dim> &p) const;

    /** get the potential value at a specified point using dealii with known cell id for the point  */
    double probe_potential(const Point<dim> &p, const int cell_index) const;

    /** Probes the field at point p that belongs in cell with cell_index. Fast, when cell_index is correct */
    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index) const;

    /** Obtain potential and electric field values in selected nodes */
    void potential_efield_at(vector<double> &potentials, vector<Tensor<1, dim>> &fields,
            const vector<int> &cells, const vector<int> &verts) const;

    /** Calculate charge densities at given nodes in given cells */
    void charge_dens_at(vector<double> &charge_dens, const vector<int> &cells, const vector<int> &verts);

    /** Run Conjugate-Gradient solver to solve matrix equation */
    int solve() { return this->solve_cg(conf->n_cg, conf->cg_tolerance, conf->ssor_param); }

    /** Setup system for solving Poisson equation */
    void setup(const double field, const double potential);

    /** Assemble the matrix equation to solve Laplace or Poisson equation
     * by appling Neumann BC (constant field) or Dirichlet BC (constant potential) on top of simubox */
    void assemble(const bool first_time);

private:
    const ParticleSpecies* particles;
    const Config::Field* conf;   ///< solver parameters
    const LinearHexahedra* interpolator;

    double applied_field;     ///< applied electric field on top of simubox
    double applied_potential; ///< applied potential on top of simubox

    typedef typename DealSolver<dim>::LinearSystem LinearSystem;
    typedef typename DealSolver<dim>::ScratchData ScratchData;
    typedef typename DealSolver<dim>::CopyData CopyData;

    double probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;
    
    double probe_efield_norm(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Write the electric potential and field to a file in vtk format */
    void write_vtk(ofstream& out) const;

    /** Mark different regions in mesh */
    void mark_mesh();

    /** Return the boundary condition value at the centroid of face */
    double get_face_bc(const unsigned int face) const;;

    /** Assemble left-hand-side of matrix equation in a parallel manner */
    void assemble_parallel();

    /** Calculate the contribution of one cell into global matrix and rhs vector */
    void assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
            ScratchData &scratch_data, CopyData &copy_data) const;

    /** Add to the right-hand-side vector for point charges, as used in PIC. */
    void assemble_space_charge();
    void assemble_space_charge_fast();
};

} // namespace femocs

#endif /* LAPLACE_H_ */
