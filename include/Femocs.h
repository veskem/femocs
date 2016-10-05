/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include "AtomReader.h"
#include "SolutionReader.h"

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
        string external_msg;  //!< path to input script
        string mesher;          //!< simulation cell finite element mesher
        string mesh_quality;    //!< the minimum quality (maximum radius-edge ratio) of tetrahedra
        string infile;          //!< path to input script
        double latconst;        //!< lattice constant
        double coord_cutoff;    //!< cutoff distance in Angstroms for Coordination analysis
        int nnn;                //!< number of nearest neighbours for given crystal structure
        int nt;                 //!< number of OpenMP threads
        double neumann;         //!< value of Neumann boundary condition
        bool postprocess_marking; //!< make extra effort to mark correctly the vacuum nodes in shadowed area
        bool refine_apex;       //!< add elements to the nanotip apex

        /**  minimum distance between atoms from current and previous run so that
         * their movement is considered to be sufficiently big to recalculate electric field */
        double significant_distance;

        //!< Radius of cylinder where surface atoms are not coarsened; zero enables coarsening of all atoms.
        double rmin_coarse;
        //!< Radius of cylinder out of which the surface atoms are coarsened with constant cutoff.
        double rmax_coarse;
        //!< Factor that is proportional to the extent of surface coarsening; zero turns coarsening off.
        double coarse_factor;
        //!< Distance from surface edge where atoms are picked for rectangularization
        double rmin_rectancularize;
        //!< Width of moving average while smoothing the electric field; 0 turns smoothing off
        double movavg_width;
    };

    Config conf;          //!< Femocs configuration parameters

    /**
     * The function to generate FEM mesh and to solve differential equation(s).
     * @param E_field - long range electric field
     */
    const void run(double E_field, string msg);

    const void import_atoms(int n_atoms, const double* coordinates, const double* box,
            const int* nborlist);
    const void import_atoms(int n_atoms, double* x, double* y, double* z, int* types);
    const void import_atoms(string file_name);
    const void export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    const void export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm,
            const int* nborlist);

private:
    bool solution_valid;
    femocs::AtomReader reader;
    femocs::SolutionReader solution;
};

const void femocs_speaker(string path);

#endif /* FEMOCS_H_ */
