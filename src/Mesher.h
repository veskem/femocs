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
#include "Femocs.h"
#include "Surface.h"

namespace femocs {
class Vacuum;
} /* namespace femocs */

using namespace std;

namespace femocs {

/**
 * Class to create and handle FEM mesh in tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/
 */
class Mesher {
public:
    Mesher(string mesher);
    virtual ~Mesher() {
    }
    ;
    void generate_mesh(const Femocs::SimuCell* cell, shared_ptr<Surface> bulk,
            shared_ptr<Surface> surf, Vacuum* vacuum, string cmd);
    void clean_mesh(string cmd, double rmax);
    void mark_mesh(const Femocs::SimuCell* cell, const string cmd);
    void write_faces(const string file_name);
    void write_elems(const string file_name);
    void refine(const string cmd);

private:
    tetgenio tetgenIn;
    tetgenio tetgenOut;
    tetgenbehavior tetgenbeh;

    void mark_faces(const Femocs::SimuCell* cell, tetgenio* tetIO);
    void mark_elems(const Femocs::SimuCell* cell, tetgenio* tetIO);

    void write_vtk(const string file_name, const int nnodes, const int ncells, const int nmarkers,
            const REAL* nodes, const int* cells, const char* markerstr, const int celltype,
            const int nnodes_per_cell);

    bool on_face(const double face1_x, const double face1_y, const double face1_z,
            const double face2_x, const double face2_y, const double face2_z);
    void add_point(const double x, const double y, const double z, const int pmarker, const int i,
            tetgenio* tetIO);
    void generate_surface_mesh(shared_ptr<Surface> surf, tetgenio* tetIO);
    void generate_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum, tetgenio* tetIO);
    void clean_elements(tetgenio* tetIO, double rmax);
    void clean_faces(tetgenio* tetIO, double rmax);
    void clean_edges(tetgenio* tetIO, double rmax);
    void calculate(string cmd, tetgenio* t1, tetgenio* t2);
    void update_list(int* list, int* list_size, vector<bool> is_quality, int M);
};

} /* namespace femocs */

#endif /* MESHER_H_ */
