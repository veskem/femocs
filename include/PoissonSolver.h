/*
 * laplace.h -> Laplace.h
 *
 *  Created on: Jul 26, 2016
 *      Author: kristjan, Mihkel
 */

#ifndef LAPLACE_H_
#define LAPLACE_H_

#include "DealSolver.h"
#include "ParticleSpecies.h"
#include "Config.h"
#include "InterpolatorCells.h"

namespace fch {

using namespace dealii;
using namespace femocs;
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

    /**
     * method to obtain the electric potential values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return potential values in the specified nodes
     */
    vector<double> get_potential(const vector<int> &cell_indexes, const vector<int> &vert_indexes);

    /**
     * method to obtain the electric field values in selected nodes
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return electric field vectors in the specified nodes
     */
    vector<Tensor<1, dim>> get_efield(const vector<int> &cell_indexes, const vector<int> &vert_indexes) const;

    /** Calculate charge densities at given nodes in given cells */
    vector<double> get_charge_dens(const vector<int> &cell_indexes, const vector<int> &vert_indexes);

    /** Run Conjugate-Gradient solver to solve matrix equation */
    unsigned int solve() { return this->solve_cg(conf->n_phi, conf->phi_error, conf->ssor_param); }

    /** Setup system for solving Poisson equation */
    void setup(const double field, const double potential=0);

    /** Assemble the matrix equation to solve Laplace equation
     * by appling Neumann BC (constant field) on top of simubox */
    void assemble_laplace(const bool first_time);

    /** Assemble the matrix equation to solve Poisson equation
     * by appling Neumann BC (constant field) or Dirichlet BC (constant potential) on top of simubox
     * as specified in config file. */
    void assemble_poisson(const bool first_time, const bool write_time);

private:
    const ParticleSpecies* particles;
    const Config::Field* conf;   ///< solver parameters
    const LinearHexahedra* interpolator;

    double applied_field;     ///< applied electric field on top of simubox
    double applied_potential; ///< applied potential on top of simubox

    double probe_potential(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;
    
    double probe_efield_norm(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    Tensor<1, dim, double> probe_efield(const Point<dim> &p, const int cell_index, Mapping<dim,dim>& mapping) const;

    /** Write the electric potential and field to a file in vtk format */
    void write_vtk(ofstream& out) const;

    /** Mark different regions in mesh */
    void mark_mesh();

    /** Return the boundary condition value at the centroid of face */
    double get_face_bc(const unsigned int face) const;

    /** @brief Reset the system and assemble the LHS matrix
     * Calculate sparse matrix elements
     * according to the Laplace equation weak formulation
     * This should be the first function call to setup the equations (after setup_system() ).
     */
    void assemble_lhs();

    /** Add to the right-hand-side vector for point charges, as used in PIC. */
    void assemble_space_charge();
    void assemble_space_charge_fast();
};

} // namespace fch

#endif /* LAPLACE_H_ */
