/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include "Macros.h"
#include "AtomReader.h"
#include "Primitives.h"

using namespace std;

/**
 * Main class to hold Femocs object
 */
class Femocs {
public:
    /**
     * Constructor of Femocs reads in and stores input parameters.
     * @param file_name - path to Femocs input script.
     */
    Femocs(string file_name);

    /** Femocs destructor */
    ~Femocs();

    /** Struct holding data about input parameters. */
    struct Config {
        string mesher;          //!< simulation cell finite element mesher
        string infile;          //!< path to input script
        double latconst;        //!< lattice constant
        double coord_cutoff;    //!< cutoff distance in Angstroms for Coordination analysis
        int nnn;                //!< number of nearest neighbours for given crystal structure
        int nt;                 //!< number of OpenMP threads
        int poly_degree;        //!< polynomial degree of the finite elements (1-linear, 2-quadratic, ...)
        double neumann;         //!< value of Neumann boundary condition
        double zmax_coarse;     //!< maximum z-coordinate of hollow cylinder where coarsening of surface atoms takes place
        double rmin_coarse;     //!< inner radius of hollow cylinder where coarsening of surface atoms takes place
        bool coarsen;           //!< flag controlling whether to coarsen or not the surface atoms
        string mesh_quality;    //!< the minimum quality (maximum radius-edge ratio) of tetrahedra
    };

    Config conf;          //!< Femocs configuration parameters

    /**
     * The function to get input data, to start the Femocs simulation and
     * to return simulation results to the caller.
     * @param E0 - long range electric field
     * @param BC - boundary conditions
     * @param phi_guess - guess values of field potential for FEM solver
     * @param grid_spacing - FDM grid spacing in x, y and z direction
     */
    const void run(double E_field, double*** phi);

    const void import_atoms(int n_atoms, double* x, double* y, double* z, int* types);

private:
    femocs::AtomReader reader;
};

const void femocs_speaker(string path);

#endif /* FEMOCS_H_ */
