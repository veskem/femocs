/*
 * Mesher.h
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#ifndef MESHER_H_
#define MESHER_H_

#include "Medium.h"
#include "AtomReader.h"
#include "../lib/tetgen.h"
#include <string>
#include <memory>
using namespace std;

namespace femocs {

/**
 * Class to create and handle FEM mesh in tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/
 */
class Mesher {
public:
	Mesher(string);
	virtual ~Mesher(){};
	void generate_mesh(shared_ptr<Surface>, string);
	void clean_mesh(string, double);
	void mark_mesh(AtomReader::Data, string);
	void smooth_mesh(int);
	void output_mesh();

private:
	tetgenio tetgenIn;
	tetgenio tetgenOut;
	tetgenbehavior tetgenbeh;
	void clean_elements(tetgenio*, double);
	void clean_faces(tetgenio*, double);
	void clean_edges(tetgenio*, double);
	void calculate(string, tetgenio*, tetgenio*);
	void update_list(int*, int*, vector<bool>, int);
};

} /* namespace femocs */

#endif /* MESHER_H_ */
