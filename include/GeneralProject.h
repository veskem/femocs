/*
 * GeneralProject.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef GENERALPROJECT_H_
#define GENERALPROJECT_H_

#include "AtomReader.h"
#include "Config.h"
#include "TetgenMesh.h"

using namespace std;
namespace femocs {

/**
 * Template for creating FEM project for FEMOCS
 * See other Project* classes to get inspiration for building up the work flow.
 */
class GeneralProject {
public:
    GeneralProject();
    GeneralProject(AtomReader &reader, Config& conf);
    virtual ~GeneralProject() {};

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @return  0 - function completed normally; 1 - function did not complete normally
     */
    virtual int run(const int timestep=-1, const double time=-1) = 0;

    /** Export the solution data in the location of imported atoms */
    virtual int export_data(double* data, const int n_points, const string& data_type) = 0;

    /** Export the integer data in the location of imported atoms */
    virtual int export_data(int* data, const int n_points, const string& data_type) = 0;

    /** Export pointer to double data and nr of data clusters in the pointer */
    virtual int export_data(const double** data, const string& data_type) const;

    /** Export pointer to integer data and nr of data clusters in the pointer */
    virtual int export_data(const int** data, const string& data_type) const;

    /** Interpolate the solution data in the location of specified points */
    virtual int interpolate(double* data, int* flag,
            const int n_points, const string& data_type, const bool near_surface,
            const double* x, const double* y, const double* z) = 0;

    /** Read and generate simulation data to continue running interrupted simulation */
    virtual int restart(const string& path_to_file) = 0;

protected:
    // references instead of pointers to make their access more convenient
    AtomReader &reader;    ///< All the imported atoms
    Config &conf;          ///< Configuration parameters

    TetgenMesh mesh1;      ///< FEM mesh for the whole simulation domain
    TetgenMesh mesh2;
    TetgenMesh *new_mesh;  ///< Pointer to mesh where the new one will be generated
    TetgenMesh *mesh;      ///< Pointer to readily available mesh

private:
    // needed to give initial values for references
    // in case of call for empty constructor
    AtomReader dummy_reader;
    Config dummy_config;
};

} /* namespace femocs */

#endif /* GENERALPROJECT_H_ */
