/*
 * GeneralProject.h
 *
 *  Created on: 8.2.2018
 *      Author: veske
 */

#ifndef GENERALPROJECT_H_
#define GENERALPROJECT_H_

namespace femocs {

#include "AtomReader.h"
#include "Config.h"
#include "Coarseners.h"
#include "SolutionReader.h"
#include "TetgenMesh.h"

/*
 *
 */
class GeneralProject {

    GeneralProject(const AtomReader &reader, const Config& conf);

    virtual ~GeneralProject() = 0;

    /** Generate meshes using the imported atomistic data */
    virtual int generate_mesh();

    /** Function to generate FEM mesh and to solve differential equation(s)
     * @return  0 - function completed normally; 1 - function did not complete normally
     */
    virtual int run();

    /** Force the data to the files for debugging purposes */
    virtual int force_output();

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
};

} /* namespace femocs */

#endif /* GENERALPROJECT_H_ */
