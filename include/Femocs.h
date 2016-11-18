/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include "AtomReader.h"
#include "Config.h"
#include "Interpolator.h"
#include "TetgenMesh.h"
#include "SolutionReader.h"
#include "Media.h"

using namespace std;
namespace femocs {

/** Main class to hold Femocs object */
class Femocs {
public:
    /**
     * Femocs constructor reads and stores configuration parameters
     * @param path_to_conf      path to input script holding the configuration parameters
     */
    Femocs(string path_to_conf);

    /** Femocs destructor */
    ~Femocs();

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param E_field   long range electric field strength
     * @param msg       message from the host: time step, file name etc
     */
    const int run(double E_field, string msg);

    /** Function to import atoms from PARCAS
     * @param n_atoms       number of imported atoms
     * @param coordinates   normalised coordinates of atoms in PARCAS format
     * @param box           size on simulation box in Angstroms
     * @param nborlist      neighbour list for atoms
     */
    const int import_atoms(int n_atoms, double* coordinates, double* box, int* nborlist);

    /** Function to import coordinates of atoms
     * @param n_atoms   number of imported atoms
     * @param x         x-coordinates of the atoms
     * @param y         y-coordinates of the atoms
     * @param z         z-coordinates of the atoms
     * @param types     types of the atoms
     */
    const int import_atoms(int n_atoms, double* x, double* y, double* z, int* types);

    /** Function to import atoms from file
     * @param file_name path to the file
     */
    const int import_atoms(const string& file_name);
    
    /** Function to export the calculated electric field
     * @param n_atoms   number of points where electric field was calculted
     * @param Ex        x-component of electric field
     * @param Ey        y-component of electric field
     * @param Ez        z-component of electric field
     * @param Enorm     norm of electric field
     */
    const int export_elfield(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    
    /** Function to linearly interpolate electric field at given points
     * @param n_points  number of points where electric field is interpolated
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     * @param flag      index of first point outside the mesh; index is 1-based; flag == 0 means all the points were in the mesh
     */
    const int interpolate_elfield(int n_points, double* x, double* y, double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Function to linearly interpolate electric potential at given points
     * @param n_points  number of points where electric potential is interpolated
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param phi       electric potential
     * @param flag      index of first point outside the mesh; index is 1-based; flag == 0 means all the points were in the mesh
     */
    const int interpolate_phi(int n_points, double* x, double* y, double* z, double* phi, int* flag);

    const int parse_command(const string & command, int* arg);

    const int parse_command(const string & command, double* arg);

private:
    string home;
    bool skip_calculations;
    AtomReader reader;
    Config conf;
    Surface dense_surf;
    TetgenMesh tetmesh_vacuum;
    TetgenMesh tetmesh_bulk;
    SolutionReader solution = SolutionReader(&tetmesh_vacuum);
    Interpolator interpolator = Interpolator(&solution);
};

} /* namespace femocs */

const void femocs_speaker(string path);

#endif /* FEMOCS_H_ */
