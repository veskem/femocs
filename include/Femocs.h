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
     * by using the parameters specified in configuration script.
     * @param timestep  active time step in the host code; if not provided, internal counter will be used
     * @return          0 - function completed normally; 1 - function did not complete normally
     */
    int run(const int timestep=-1, const double time=-1);

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
     * @param data       array where solution data is written; vector data is written component-wise, i.e in a from x1,y1,z1,x2,y2...
     * @param n_points   number of first imported points where solution is exported; <= 0 turns export off
     * @param data_type  label of data to be exported
     */
    int export_data(double* data, const int n_points, const string& data_type);

    /** Interpolate the solution data in the location of specified points
     * @param data          array where solution data is written; vector data is written component-wise, i.e in a from x1,y1,z1,x2,y2...
     * @param flag          indicators showing the location of point; 0 - point was inside the mesh, 1 - point was outside the mesh
     * @param n_points      number of points in x,y,z arrays; <= 0 turns interpolation off
     * @param data_type     label of data to be exported
     * @param near_surface  data points are located near the surface
     * @param x,y,z         coordinates of the points
     */
    int interpolate(double* data, int* flag,
            const int n_points, const string& data_type, const bool near_surface,
            const double* x, const double* y, const double* z);

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

private:
    double t0;

    Config conf;             ///< configuration parameters
    AtomReader reader;       ///< all the imported atoms
    GeneralProject *project; ///< project Femocs is going to run

    void perform_full_analysis(const int* nborlist);
    void perform_pseudo_analysis();
};

} /* namespace femocs */

#endif /* FEMOCS_H_ */
