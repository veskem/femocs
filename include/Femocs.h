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
#include "TetgenMesh.h"
#include "SolutionReader.h"

using namespace std;
namespace femocs {

/** Main class to hold Femocs object */
class Femocs {
public:
    /**
     * Constructor of Femocs reads in and stores input parameters
     * @param file_name - path to Femocs input script
     */
    Femocs(string file_name);

    /** Femocs destructor */
    ~Femocs();

    /** Struct holding data about input parameters. */
    struct Config {
        string mesh_quality;      ///< Minimum quality (maximum radius-edge ratio) of tetrahedra
        string infile;            ///< Path to input script
        double latconst;          ///< Lattice constant
        double coord_cutoff;      ///< Cut-off distance in Angstroms for Coordination analysis
        int nnn;                  ///< Number of nearest neighbours for given crystal structure
        int nt;                   ///< Number of OpenMP threads
        double neumann;           ///< Value of Neumann boundary condition
        bool postprocess_marking; ///< Make extra effort to mark correctly the vacuum nodes in shadowed area
        bool refine_apex;         ///< Add elements to the nanotip apex
        double zbox_above;        ///< Space added above the maximum z-coordinate of surface
        double zbox_below;        ///< Space added below the minimum z-coordinate of surface

        /// Distance from surface edge where atoms are picked for rectangularization
        double rmin_rectancularize;

        /** Minimum distance between atoms from current and previous run so that their
         * movement is considered to be sufficiently big to recalculate electric field;
         * 0 turns the check off */
        double significant_distance;

        /// Radius of cylinder where surface atoms are not coarsened; 0 enables coarsening of all atoms
        double radius;
        
        /// Factor that is proportional to the extent of surface coarsening; 0 turns coarsening off
        double coarse_factor;
        
        /// Factor that is proportional to the extent of surface smoothing; 0 turns smoothing off
        double smooth_factor;
        
        /// Width of moving average while smoothing the electric field; 0 turns smoothing off
        double movavg_width;

        /// Number of bins in smoother histogram; 1 or less turns off the histogram smoother
        int n_bins;

    } conf;

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param E_field   long range electric field strength
     * @param msg       message from the host: time step, file name etc
     */
    const void run(double E_field, string msg);

    /** Function to import atoms from PARCAS
     * @param n_atoms       number of imported atoms
     * @param coordinates   normalised coordinates of atoms in PARCAS format
     * @param box           size on simulation box in Angstroms
     * @param nborlist      neighbour list for atoms
     */
    const void import_atoms(int n_atoms, const double* coordinates, const double* box, const int* nborlist);

    /** Function to import coordinates of atoms
     * @param n_atoms   number of imported atoms
     * @param x         x-coordinates of the atoms
     * @param y         y-coordinates of the atoms
     * @param z         z-coordinates of the atoms
     * @param types     types of the atoms
     */
    const void import_atoms(int n_atoms, double* x, double* y, double* z, int* types);

    /** Function to import atoms from file
     * @param file_name path to the file
     */
    const void import_atoms(string file_name);
    
    /** Function to export the calculated electric field
     * @param n_atoms   number of points where electric field was calculted
     * @param Ex        x-component of electric field
     * @param Ey        y-component of electric field
     * @param Ez        z-component of electric field
     * @param Enorm     norm of electric field
     */
    const void export_elfield(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    
    /** Function to linearly interpolate electric field at given points
     * @param n_points  number of points where electric field is interpolated
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     */
    const void interpolate_elfield(int n_points, double* x, double* y, double* z, double* Ex,
            double* Ey, double* Ez, double* Enorm);

    /** Function to linearly interpolate electric potential at given points
     * @param n_points  number of points where electric potential is interpolated
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param phi       electric potential
     */
    const void interpolate_phi(int n_points, double* x, double* y, double* z, double* phi);

private:
    string home;
    bool solution_valid;
    AtomReader reader;
    Interpolator interpolation;
    TetgenMesh tetmesh_vacuum;
    TetgenMesh tetmesh_bulk;
    SolutionReader solution = SolutionReader(&tetmesh_vacuum);
};

} /* namespace femocs */

const void femocs_speaker(string path);

#endif /* FEMOCS_H_ */
