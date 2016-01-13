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
	~Femocs() {};

	/** Struct holding data about input parameters. */
	struct Config {
		string smoother;        //!< surface mesh smoother
		string extracter;       //!< surface extracter from bulk material
		string mesher;          //!< simulation cell finite element mesher
		double coord_cutoff;    //!< cutoff distance in Angstroms for Coordination analysis
		double tetgen_cutoff;   //!< cutoff distance in Angstroms for removing too big mesh elements
		int nnn;                //!< number of nearest neighbours for given crystal structure
		int nt;                 //!< number of OpenMP threads
	};
	const Config conf;          //!< Femocs configuration parameters

private:
	/**
	 * Function to get configuration parameters from input script
	 * @param file_name - path to input script
	 * @return configuration parameters to Femocs::Config struct
	 */
	const Config parse_input_script(const string& file_name) const;
};

} /* namespace femocs */
#endif /* FEMOCS_H_ */
