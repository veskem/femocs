/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include "AtomReader.h"
#include "Interpolator.h"

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
    Femocs(string file_name);

    /** Femocs destructor */
    ~Femocs();

    /** Struct holding data about input parameters. */
    struct Config {
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

        /** Minimum distance between atoms from current and previous run so that
         * their movement is considered to be sufficiently big to recalculate electric field;
         * zero turns the check off */
        double significant_distance;

        /** Radius of cylinder where surface atoms are not coarsened;
         * zero enables coarsening of all atoms. */
        double radius;
        
        /** Factor that is proportional to the extent of surface coarsening; 
         * zero turns coarsening off. */
        double coarse_factor;
        
        /** Factor that is proportional to the extent of surface smoothing; 
         * zero turns smoothing off. */
        double smooth_factor;
        
        /** Distance from surface edge where atoms are picked for rectangularization */
        double rmin_rectancularize;
        
        /** Width of moving average while smoothing the electric field; 
         * 0 turns smoothing off */
        double movavg_width;

        /** number of bins in smoother histogram; 
         * 1 or less turns off the histogram smoother */
        int n_bins;

        /** Space added above the maximum z-coordinate of surface */
        double zbox_above;

        /** Space added below the minimum z-coordinate of surface */
        double zbox_below;
    } conf;

    /**
     * The function to generate FEM mesh and to solve differential equation(s).
     * @param E_field - long range electric field
     */
    const void run(double E_field, string msg);

    const void import_atoms(int n_atoms, const double* coordinates, const double* box, const int* nborlist);
    const void import_atoms(int n_atoms, double* x, double* y, double* z, int* types);
    const void import_atoms(string file_name);
    
    const void export_solution(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    
    const void interpolate_solution(int n_atoms, double* x, double* y, double* z, 
                                    double* Ex, double* Ey, double* Ez, double* Enorm);

private:
    string home;
    bool solution_valid;
    AtomReader reader;
    Interpolator interpolation;
};

} /* namespace femocs */

const void femocs_speaker(string path);

#endif /* FEMOCS_H_ */
