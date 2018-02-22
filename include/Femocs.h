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
#include "GeneralProject.h"

using namespace std;
namespace femocs {

/** Main class to hold Femocs object */
class Femocs {
public:
    /**
     * Femocs constructor reads and stores configuration parameters and initialises other data
     * @param path_to_conf      path to the file holding the configuration parameters
     */
    Femocs(const string &path_to_conf);

    /** Femocs destructor */
    ~Femocs();

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @param elfield   long range electric field strength
     * @param timestep  active time step in the host code
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const double elfield, const string& timestep);

    /** Function to generate FEM mesh and to solve differential equation(s)
     * by using the parameters specified in configuration script.
     * @param timestep  active time step in the host code; if not provided, internal counter will be used
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const int timestep=-1);

    /** Function to generate artificial nanotip without crystallographic features
     * @param height      height of the tip sides in units of radius; the total tip height is therefore (height + 1)*radius; negative value makes open void instead of tip
     * @param radius      radius of the tip; if not specified, its value is taken from the configuration script
     * @param resolution  distance between the atoms; if not specified, its value will be equal to lattice constant
     * @return            success of the operation (always 0)
     */
    int generate_nanotip(const double height, const double radius=-1, const double resolution=-1);

    /** Function to import atoms from PARCAS
     * @param n_atoms       number of imported atoms
     * @param coordinates   normalised coordinates of atoms in PARCAS format
     * @param box           size on simulation box in Angstroms
     * @param nborlist      neighbour list for atoms
     * @return              0 - import succeeded and current run differs from previous one, whole field calculation will be performed;
     *                      1 - import failed or current run is similar to the previous one, field calculation will be skipped
     */
    int import_atoms(const int n_atoms, const double* coordinates, const double* box, const int* nborlist);

    /** Function to import coordinates and types of atoms
     * @param n_atoms   number of imported atoms
     * @param x         x-coordinates of the atoms
     * @param y         y-coordinates of the atoms
     * @param z         z-coordinates of the atoms
     * @param types     types of the atoms
     * @return          0 - import succeeded and current run differs from previous one, whole field calculation will be performed;
     *                  1 - import failed or current run is similar to the previous one, field calculation will be skipped
     */
    int import_atoms(const int n_atoms, const double* x, const double* y, const double* z, const int* types);

    /** Function to import atoms from file
     * @param file_name path to the file
     * @param add_noise add random noise to the imported atom coordinates to emulate real simulation
     * @return          0 - import succeeded, 1 - import failed
     */
    int import_atoms(const string& file_name, const int add_noise=0);
    
    /** Export the types of all the atoms as seen by FEMOCS
     * @param n_atoms   number of atoms to export; n_atoms <= 0 turns the export off
     * @param types     array where the atom types are written
     * @return          boolean whether there are any clustered or evaporated atom
     */
    int export_atom_types(const int n_atoms, int* types);

    /** Function to export the calculated electric field on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param Ex        x-component of electric field
     * @param Ey        y-component of electric field
     * @param Ez        z-component of electric field
     * @param Enorm     norm of electric field
     * @return          success of the operation (always 0)
     */
    int export_elfield(const int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    
    /** Function to export the calculated temperatures on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param T         temperature in the atom location
     * @return          success of the operation (always 0)
     */
    int export_temperature(const int n_atoms, double* T);

    /** Calculate and export charges & forces on imported atom coordinates
     * @param n_atoms   number of points of interest; n_atoms <= 0 turns the export off
     * @param xq        charges and forces in PARCAS format (xq[0] = q1, xq[1] = Fx1, xq[2] = Fy1, xq[3] = Fz1, xq[4] = q2, xq[5] = Fx2 etc)
     * @return          success of the operation (always 0)
     */
    int export_charge_and_force(const int n_atoms, double* xq);

    /** Export Laplace + Coulomb force and pair potential on imported atoms
     * @param n_atoms  number of first atoms the data will be exported
     * @param xnp      forces in PARCAS format & units (xnp[0] = Fx1, xnp[1] = Fy1, xnp[2] = Fz1, xnp[3] = Fx2 etc)
     * @param Epair    potential energy per atom
     * @param Vpair    total potential energy of atoms. Pot. due to Coloumb forces are added here. NOTE: Lorentz is missing!
     * @return         success of the operation (always 0)
     */
    int export_force_and_pairpot(const int n_atoms, double* xnp, double* Epair, double* Vpair);

    /** Function to linearly interpolate electric field at points anywhere in space
     * @param n_points  number of points where electric field is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Function to linearly interpolate electric field at points on or near the surface.
     * It is faster than interpolate_elfield but is inadequate for points far away from surface
     * @param n_points  number of points where electric field is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param Ex        x-component of the interpolated electric field
     * @param Ey        y-component of the interpolated electric field
     * @param Ez        z-component of the interpolated electric field
     * @param Enorm     norm of the interpolated electric field
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_surface_elfield(const int n_points, const double* x, const double* y, const double* z,
            double* Ex, double* Ey, double* Ez, double* Enorm, int* flag);

    /** Function to linearly interpolate electric potential at given points
     * @param n_points  number of points where electric potential is interpolated; n_points <= 0 turns the interpolation off
     * @param x         x-coordinates of the points of interest
     * @param y         y-coordinates of the points of interest
     * @param z         z-coordinates of the points of interest
     * @param phi       electric potential
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @return          0 - function used solution from current run; 1 - function used solution from previous run
     */
    int interpolate_phi(const int n_points, const double* x, const double* y, const double* z, double* phi, int* flag);

    /** Export the solution data in the location of imported atoms
     * @param n_points  number of first imported points where solution is exported; <= 0 turns export off
     * @param cmd       label specifying the data to be exported
     * @param data      array where solution data is written; vector data is written component-wise, i.e in a from x1,y1,z1,x2,y2...
     */
    int export_results(const int n_points, const char cmd, double* data);

    /** Export the solution data in the location of specified points
     * @param n_points  number of points in x,y,z arrays; <= 0 turns interpolation off
     * @param cmd       label specifying the data to be exported
     * @param x,y,z     coordinates of the points
     * @param data      array where solution data is written; vector data is written component-wise, i.e in a from x1,y1,z1,x2,y2...
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     */
    int interpolate_results(const int n_points, const char cmd,
            const double* x, const double* y, const double* z, double* data, int* flag);

    /** Export the solution data in the location of specified points.
     * Points are assumed to be near the surface, allowing to make interpolation faster.
     * @param n_points  number of points in x,y,z arrays; <= 0 turns interpolation off
     * @param cmd       label specifying the data to be exported
     * @param x,y,z     coordinates of the points
     * @param data      array where solution data is written; vector data is written component-wise, i.e in a from x1,y1,z1,x2,y2...
     * @param flag      indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     */
    int interpolate_surface_results(const int n_points, const char cmd,
            const double* x, const double* y, const double* z, double* data, int* flag);

    /**
     * Function to parse integer argument of the command from input script
     * @param command   name of the command which's argument should be parsed
     * @param arg       parsed argument
     * @return          0 - command was found and arg was modified; 1 - command was not found and arg was not modified
     */
    int parse_command(const string& command, int* arg);

    /**
     * Function to parse double argument of the command from input script
     * @param command   name of the command which's argument should be parsed
     * @param arg       parsed argument
     * @return          0 - command was found and arg was modified; 1 - command was not found and arg was not modified
     */
    int parse_command(const string& command, double* arg);

    /**
     * Function to parse char array argument of the command from input script
     * @param command   name of the command which's argument should be parsed
     * @param arg       parsed argument
     * @return          0 - command was found and arg was modified; 1 - command was not found and arg was not modified
     */
    int parse_command(const string& command, char* arg);

    /**
     * Function to parse string argument of the command from input script
     * @param command   name of the command which's argument should be parsed
     * @param arg       parsed argument
     * @return          0 - command was found and arg was modified; 1 - command was not found and arg was not modified
     */
    int parse_command(const string& command, string& arg);

    /** Force the data to the files for debugging purposes */
    int force_output();

private:
    double t0;

    Config conf;             ///< configuration parameters
    AtomReader reader;       ///< all the imported atoms
    GeneralProject *project; ///< project Femocs is going to run
};

} /* namespace femocs */

#endif /* FEMOCS_H_ */
