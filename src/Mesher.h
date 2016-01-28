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
    void clean_mesh(string cmd, double rmax, const Femocs::SimuCell* cell);
    void mark_mesh(const Femocs::SimuCell* cell, const string cmd);
    void add_vacuum(const string cmd, Femocs::SimuCell* cell);
    void write_faces(const string file_name);
    void write_elems(const string file_name);
    void refine(const string cmd);
    void generate_surface_mesh(shared_ptr<Surface> surf);
    void generate_bulk_mesh(const Femocs::SimuCell* cell, shared_ptr<Surface> bulk);
    void generate_volume_mesh(shared_ptr<Surface> bulk, Vacuum* vacuum);
    void unite(const string cmd, const Mesher* mesh2);
    void unite2(const string cmd, const Mesher* mesh2, const Femocs::SimuCell* cell);
    void unite3(const tetgenio* tetIO_bulk, const tetgenio* tetIO_surf, const Femocs::SimuCell* cell, const string cmd);

    void calc(const string cmd);
    tetgenio tetgenOut;

private:
    tetgenbehavior tetgenbeh;

    void mark_faces(const Femocs::SimuCell* cell, tetgenio* tetIO);
    void mark_elems(const Femocs::SimuCell* cell, tetgenio* tetIO);

    void write_vtk(const string file_name, const int nnodes, const int ncells, const int nmarkers,
            const REAL* nodes, const int* cells, const char* markerstr, const int celltype,
            const int nnodes_per_cell);

    bool on_face(const double f1n1, const double f1n2, const double f1n3, const double f2);
    void add_point(const double x, const double y, const double z, const int pmarker, const int i,
            tetgenio* tetIO);

    void clean_elements(tetgenio* tetIO, double rmax);
    void clean_faces(tetgenio* tetIO, double rmax);
    void clean_faces2(tetgenio* tetIO, double rmax, const Femocs::SimuCell* cell);
    void clean_edges(tetgenio* tetIO, double rmax);
    void calculate(string cmd, tetgenio* t1, tetgenio* t2);
    void update_list(int* list, int* list_size, vector<bool> is_quality, int M);
};

} /* namespace femocs */

#endif /* MESHER_H_ */
