/*
 * Femocs.h
 *
 *  Created on: 14.12.2015
 *      Author: veske
 */

#ifndef FEMOCS_H_
#define FEMOCS_H_

#include <string>

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
    Femocs(const string& file_name);
    ~Femocs() {
    }
    ;

    /** Struct holding data about input parameters. */
    struct Config {
        string extracter;       //!< surface extracter from bulk material
        string mesher;          //!< simulation cell finite element mesher
        double latconst;        //!< lattice constant
        double coord_cutoff;    //!< cutoff distance in Angstroms for Coordination analysis
        double tetgen_cutoff;   //!< cutoff distance in Angstroms for removing too big mesh elements
        int nnn;                //!< number of nearest neighbours for given crystal structure
        int nt;                 //!< number of OpenMP threads
        double zmaxbox;         //!< maximum z-coordinate of simulation box
    };

    /** Struct holding data about general simulation cell parameters. */
    struct SimuCell {
        double xmin;    //!< minimum x-coordinate of atoms
        double xmax;    //!< maximum x-coordinate of atoms
        double ymin;    //!< minimum y-coordinate of atoms
        double ymax;    //!< maximum y-coordinate of atoms
        double zmin;    //!< minimum z-coordinate of atoms
        double zmax;    //!< maximum z-coordinate of atoms
        double zmaxbox; //!< maximum z-coordinate of simulation box
        int type_bulk = 1; //!< type of bulk material
        int type_surf = 2; //!< type of open material surface
        int type_vacancy = 3; //!< type of vacancies
        int type_vacuum = 3;  //!< type of vacuum
        int type_fixed = -1;  //!< type of fixed atoms
    };

    Config conf;          //!< Femocs configuration parameters
    SimuCell simucell;    //!< General data about simulation cell

private:
    /**
     * Function to get configuration parameters from input script
     * @param file_name - path to input script
     * @return configuration parameters to Femocs::Config struct
     */
    const Config parse_input_script(const string& file_name) const;

    /** Function to initialise SimuCell values */
    const SimuCell init_simucell() const;
};

} /* namespace femocs */
#endif /* FEMOCS_H_ */
