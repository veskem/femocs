/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/base/timer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <memory>
#include <string>

#define DEBUGMODE true
#define LIBARYMODE false
#define HELMODMODE false
#define KIMOCSMODE false

using namespace std;
namespace femocs {

/**
 * Main class to hold Femocs object
 */
class Femocs {
public:
    /**
     * Constructor of Femocs reads in and stores input parameters.
     * @param file_name - path to Femocs input script.
     */
    Femocs(const string& file_name);

    /** Struct holding data about input parameters. */
    struct Config {
        string extracter;       //!< surface extracter from bulk material
        string mesher;          //!< simulation cell finite element mesher
        string infile;          //!< path to input script
        double latconst;        //!< lattice constant
        double coord_cutoff;    //!< cutoff distance in Angstroms for Coordination analysis
        double tetgen_cutoff;   //!< cutoff distance in Angstroms for removing too big mesh elements
        int nnn;                //!< number of nearest neighbours for given crystal structure
        int nt;                 //!< number of OpenMP threads
        double zmaxbox;         //!< maximum z-coordinate of simulation box
        int poly_degree;        //!< polynomial degree of the finite element (1-linear, 2-quadratic,..)
        double neumann;         //!< value of Neumann boundary condition
    };

    /** Struct holding data about general simulation cell parameters. */
    struct SimuCell {
        double xmin;    //!< minimum x-coordinate of atoms
        double xmax;    //!< maximum x-coordinate of atoms
        double ymin;    //!< minimum y-coordinate of atoms
        double ymax;    //!< maximum y-coordinate of atoms
        double zmin;    //!< minimum z-coordinate of atoms
        double zmax;    //!< maximum z-coordinate of atoms
        double zminbox; //!< minimum z-coordinate of simulation box
        double zmaxbox; //!< maximum z-coordinate of simulation box
        int type_bulk = 1; //!< type of bulk material
        int type_surf = 2; //!< type of open material surface
        int type_vacancy = 3; //!< type of vacancies
        int type_vacuum = 3;  //!< type of vacuum
        int type_edge = 0;  //!< type of the rim/outer edge of surface
        int type_fixed = -1;  //!< type of fixed atoms
        int type_xmin = 4;  //!< type of atom on negative x-face of simulation cell
        int type_ymin = 5;  //!< type of atom on negative y-face of simulation cell
        int type_zmin = 6;  //!< type of atom on negative z-face of simulation cell
        int type_xmax = 10;  //!< type of atom on positive x-face of simulation cell
        int type_ymax = 9;  //!< type of atom on positive y-face of simulation cell
        int type_zmax = 8;  //!< type of atom on positive z-face of simulation cell
        int type_none = 7;  //!< type of atom with unknown position
    };

    Config conf;          //!< Femocs configuration parameters
    SimuCell simucell;    //!< General data about simulation cell

    /**
     * The function to get input data, to start the Femocs simulation and
     * to return simulation results to the caller.
     * @param E0 - long range electric field
     * @param BC - boundary conditions
     * @param phi_guess - guess values of field potential for FEM solver
     * @param grid_spacing - FDM grid spacing in x, y and z direction
     */
    const void run_femocs(const double E0, double*** BC, double*** phi_guess, const double* grid_spacing);

private:
    /**
     * Function to get configuration parameters from input script
     * @param file_name - path to input script
     * @return configuration parameters to Femocs::Config struct
     */
    const Config parse_input_script(const string& file_name) const;

    /** Function to initialise SimuCell values */
    const SimuCell init_simucell() const;

    /** Function to run the simulation commands */

};

} /* namespace femocs */
#endif /* FEMOCS_H_ */
