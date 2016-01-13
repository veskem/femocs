/*
 * Mesher.h
 *
 *  Created on: 16.12.2015
 *      Author: veske
 */

#ifndef MESHER_H_
#define MESHER_H_

#include <memory>
#include <string>
#include <vector>

#include "../lib/tetgen.h"
#include "AtomReader.h"
#include "Medium.h"

using namespace std;

namespace femocs {

/**
 * Class to create and handle FEM mesh in tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/
 */
class Mesher {
public:
    Mesher(string mesher);
    virtual ~Mesher() {};
    void generate_mesh(shared_ptr<Surface> surf, string cmd);
    void clean_mesh(string cmd, double rmax);
    void mark_mesh(AtomReader::Data data, string cmd);
    void smooth_mesh(int nr_of_iterations);
    void output_mesh();

private:
    tetgenio tetgenIn;
    tetgenio tetgenOut;
    tetgenbehavior tetgenbeh;
    void clean_elements(tetgenio* tetIO, double rmax);
    void clean_faces(tetgenio* tetIO, double rmax);
    void clean_edges(tetgenio* tetIO, double rmax);
    void calculate(string cmd, tetgenio* t1, tetgenio* t2);
    void update_list(int* list, int* list_size, vector<bool> is_quality, int M);
};

} /* namespace femocs */

#endif /* MESHER_H_ */
