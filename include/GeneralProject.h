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
#include "Coarseners.h"
#include "TetgenMesh.h"
#include "SolutionReader.h"

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
    virtual int run(const int timestep=-1) = 0;

    /** Force the data to the files for debugging purposes */
    virtual int force_output() = 0;

    virtual int export_results(const int n_points, const string &data_type, double* data) = 0;

    virtual int interpolate_results(const int n_points, const string &data_type, const bool surface,
            const double* x, const double* y, const double* z, double* data, int* flag) = 0;

    FieldReader fields;       ///< fields & potentials on surface atoms
    HeatReader  temperatures; ///< temperatures & current densities on bulk atoms
    ForceReader forces;       ///< forces & charges on surface atoms

protected:
    AtomReader &reader;     ///< all the imported atoms
    Config &conf;           ///< configuration parameters

    Coarseners coarseners;  ///< atomistic coarsening data & routines

    TetgenMesh mesh1;       ///< FEM mesh for the whole simulation domain
    TetgenMesh mesh2;
    TetgenMesh *new_mesh;   ///< Pointer to mesh where the new one will be generated
    TetgenMesh *mesh;       ///< Readily available mesh

private:
    AtomReader dummy_reader;
    Config dummy_config;
};

} /* namespace femocs */

#endif /* GENERALPROJECT_H_ */
