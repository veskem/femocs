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
    EmissionSolver(Triangulation<dim> *tria);
    virtual ~EmissionSolver() {}

    /** Solve the matrix equation using conjugate gradient method */
    int solve() { return this->solve_cg(conf->n_cg, conf->cg_tolerance, conf->ssor_param); }

    /** Set the boundary condition values on dofs and on centroids of surface faces.
     * The values of prev_ch_values must be the same as for dofs.
     * The values of bc_values must be on the centroids of the vacuum-material boundary faces
     * in the order specified in the get_surface_nodes() method.
     */
    void set_bcs(vector<double>* prev_ch_values, vector<double>* bc_values) {
        this->prev_ch_values = prev_ch_values;
        this->bc_values = bc_values;
    }

    /** Set the pointers for obtaining external data */
    void set_dependencies(PhysicalQuantities *pq, const Config::Heating *conf) {
        this->pq = pq;
        this->conf = conf;
    }

    /** Return the boundary condition value at the centroid of face */
    double get_face_bc(const unsigned int face) const {
        require(face < bc_values->size(), "Invalid index: " + d2s(face));
        return (*bc_values)[face];
    }

protected:
    PhysicalQuantities *pq;         ///< object to evaluate tabulated physical quantities (sigma, kappa, gtf emission)
    const Config::Heating *conf;    ///< solver parameters

    vector<double>* prev_ch_values; ///< previous current/heat values for heat/current solver
    vector<double>* bc_values;      ///< current/heat values on the centroids of surface faces for current/heat solver
    
    friend class CurrentHeatSolver<dim> ;
};

template<int dim>
class CurrentSolver : public EmissionSolver<dim> {
public:
    CurrentSolver();
    CurrentSolver(Triangulation<dim> *tria, const HeatSolver<dim> *hs);

    /** @brief assemble the matrix equation for current density calculation.
     * Calculate sparse matrix elements and right-hand-side vector
     * according to the continuity equation weak formulation and to the boundary conditions.
     */
    void assemble();
    
private:
    const HeatSolver<dim>* heat_solver;
    
    // TODO figure out what is written
    void write_vtk(ofstream& out) const;

    /** Assemble left-hand-side of matrix equation */
    void assemble_lhs();

    friend class CurrentHeatSolver<dim> ;
};

template<int dim>
class HeatSolver : public EmissionSolver<dim> {
public:
    HeatSolver();
    HeatSolver(Triangulation<dim> *tria, const CurrentSolver<dim> *cs);

    /** Assemble the matrix equation for temperature calculation
     * using Crank-Nicolson or implicit Euler time integration method. */
    void assemble(const double delta_time);

private:
    static constexpr double cu_rho_cp = 3.4496e-24;      ///< volumetric heat capacity of copper [J/(K*Ang^3)]
    const CurrentSolver<dim>* current_solver;

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

    /** Output the temperature [K] and electrical conductivity [1/(Ohm*nm)] in vtk format */
    void write_vtk(ofstream& out) const;

    friend class CurrentHeatSolver<dim> ;
};

template<int dim>
class CurrentHeatSolver : public DealSolver<dim> {
public:
    CurrentHeatSolver();
    CurrentHeatSolver(PhysicalQuantities *pq_, const Config::Heating *conf_);

    /**
     * Method to obtain the temperature values in selected nodes.
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return temperature  ///< default temperature applied on bottom of the materiale values in the specified nodes
     */
    vector<double> get_temperature(const vector<int> &cell_indexes, const vector<int> &vert_indexes);

    /**
     * Method to obtain the current density values in selected nodes.
     * @param cell_indexes global cell indexes, where the corresponding nodes are situated
     * @param vert_indexes the vertex indexes of the nodes inside the cell
     * @return current density vectors in the specified nodes
     */
    vector<Tensor<1,dim>> get_current(const vector<int> &cell_indexes, const vector<int> &vert_indexes);

    /** Set the pointers for obtaining external data */
    void set_dependencies(PhysicalQuantities *pq_, const Config::Heating *conf_);

    void setup(const double temperature);

    HeatSolver<dim> heat;        ///< data and operations of heat equation solver
    CurrentSolver<dim> current;  ///< data and operations for finding current density in material

private:
    PhysicalQuantities *pq;
    const Config::Heating *conf;

    /** Mark different regions of the mesh */
    void mark_mesh();
};

} // namespace femocs

#endif /* CURRENTSANDHEATING_H_ */
