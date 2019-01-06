/*
 * currents_and_heating.h -> CurrentsAndHeating.h
 *
 *  Created on: Jul 28, 2016
 *      Author: kristjan, Mihkel
 */

#ifndef CURRENTSANDHEATING_H_
#define CURRENTSANDHEATING_H_


#include "DealSolver.h"
#include "PhysicalQuantities.h"
#include "Config.h"
#include "PoissonSolver.h"


namespace femocs {

// forward declare some classes to make them available for declaring dependencies
template<int dim> class HeatSolver;
template<int dim> class CurrentHeatSolver;

using namespace dealii;
using namespace std;

template<int dim>
class EmissionSolver : public DealSolver<dim> {
public:
    EmissionSolver();
    EmissionSolver(Triangulation<dim> *tria, vector<double>* bc_values);
    virtual ~EmissionSolver() {}

    /** Solve the matrix equation using conjugate gradient method */
    int solve() { return this->solve_cg(conf->n_cg, conf->cg_tolerance, conf->ssor_param); }

    /** Set the pointers for obtaining external data */
    void set_dependencies(PhysicalQuantities *pq, const Config::Heating *conf) {
        this->pq = pq;
        this->conf = conf;
    }

    /** Set the boundary condition values on centroids of surface faces.
     * The values of bc_values must be on the centroids of the vacuum-material boundary faces
     * in the order specified in the get_surface_nodes() method.
     */
    void set_bcs(vector<double>* bc_values) {
        this->bc_values = bc_values;
    }

    /** Return the boundary condition value at the centroid of face */
    double get_face_bc(const unsigned int face) const {
        require(face < bc_values->size(), "Invalid index: " + d2s(face));
        return (*bc_values)[face];
    }

protected:
    PhysicalQuantities *pq;         ///< object to evaluate tabulated physical quantities (sigma, kappa, gtf emission)
    const Config::Heating *conf;    ///< solver parameters

    vector<double>* bc_values;      ///< current/heat values on the centroids of surface faces for current/heat solver

    friend class CurrentHeatSolver<dim> ;
};

template<int dim>
class CurrentSolver : public EmissionSolver<dim> {
public:
    CurrentSolver();
    CurrentSolver(Triangulation<dim> *tria, const HeatSolver<dim> *hs, vector<double>* bc_values);

    /** @brief assemble the matrix equation for current density calculation.
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the continuity equation weak formulation and to the boundary conditions.
     */
    void assemble();

private:
    const HeatSolver<dim>* heat_solver;
    
    typedef typename DealSolver<dim>::LinearSystem LinearSystem;
    typedef typename DealSolver<dim>::ScratchData ScratchData;
    typedef typename DealSolver<dim>::CopyData CopyData;

    // TODO figure out what is written
    void write_vtk(ofstream& out) const;

    /** Calculate the contribution of one cell into global matrix and rhs vector */
    void assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
            ScratchData &scratch_data, CopyData &copy_data) const;

    friend class CurrentHeatSolver<dim> ;
};

template<int dim>
class HeatSolver : public EmissionSolver<dim> {
public:
    HeatSolver();
    HeatSolver(Triangulation<dim> *tria, const CurrentSolver<dim> *cs, vector<double>* bc_values);

    /** Assemble the matrix equation for temperature calculation
     * using Crank-Nicolson or implicit Euler time integration method. */
    void assemble(const double delta_time);

private:
    // TODO shouldn't it be temperature dependent?
    static constexpr double cu_rho_cp = 3.4496e-24;  ///< volumetric heat capacity of copper [J/(K*Ang^3)]
    double one_over_delta_time;      ///< inverse of heat solver time step [1/sec]
    Vector<double> joule_heat;       ///< integral Joule heat at dofs [Watt]
    const CurrentSolver<dim>* current_solver;

    typedef typename DealSolver<dim>::LinearSystem LinearSystem;
    typedef typename DealSolver<dim>::ScratchData ScratchData;
    typedef typename DealSolver<dim>::CopyData CopyData;

    /** @brief assemble the matrix equation for temperature calculation using Crank-Nicolson time integration method
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the time dependent heat equation weak formulation and to the boundary conditions.
     */
    void assemble_crank_nicolson(const double delta_time);

    /** @brief assemble the matrix equation for temperature calculation using implicit Euler time integration method
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the time dependent heat equation weak formulation and to the boundary conditions.
     */
    void assemble_euler_implicit(const double delta_time);

    /** Calculate the contribution of one cell into global matrix and rhs vector */
    void assemble_local_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
            ScratchData &scratch_data, CopyData &copy_data) const;

    /** Output the temperature [K] and electrical conductivity [1/(Ohm*nm)] in vtk format */
    void write_vtk(ofstream& out) const;

    friend class CurrentHeatSolver<dim> ;
};

// Forward declaration to avoid including header
class EmissionReader;

template<int dim>
class CurrentHeatSolver : public DealSolver<dim> {
public:
    CurrentHeatSolver();
    CurrentHeatSolver(PhysicalQuantities *pq_, const Config::Heating *conf_, EmissionReader *emission);

    /** Obtain the temperature, potential and current density values in selected nodes. */
    void temp_phi_rho_at(vector<double> &temp, vector<double> &phi, vector<Tensor<1,dim>> &rho,
            const vector<int> &cells, const vector<int> &verts) const;

    vector<Tensor<1, dim>> get_temp_grad(const vector<int> &cell_indexes, const vector<int> &vert_indexes);

    /** Set the pointers for obtaining external data */
    void set_dependencies(PhysicalQuantities *pq_, const Config::Heating *conf_);

    void setup(const double temperature);

    HeatSolver<dim> heat;        ///< data and operations of heat equation solver
    CurrentSolver<dim> current;  ///< data and operations for finding current density in material

private:
    PhysicalQuantities *pq;
    const Config::Heating *conf;

    /** Specify allowed types of writable files */
    bool valid_extension(const string &ext) const {
        return ext == "xyz" || ext == "movie";
    }

    /** Write nodal data as xyz file */
    void write_xyz(ofstream& out) const;

    /** Mark different regions of the mesh */
    void mark_mesh();
};

} // namespace femocs

#endif /* CURRENTSANDHEATING_H_ */
