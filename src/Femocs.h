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
 * Main class of Femocs object
 */
class Femocs {
public:
	Femocs(const string&);
	~Femocs() {};
	struct Config {
		string smoother, extracter, mesher;
		double coord_cutoff, tetgen_cutoff;
		int nnn, nt;
	};
	const Config conf;

private:
	const Config parseInputFile(const string&) const;
};

} /* namespace femocs */
#endif /* FEMOCS_H_ */
