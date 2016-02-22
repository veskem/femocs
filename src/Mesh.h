/*
 * Mesh.h
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#ifndef MESH_H_
#define MESH_H_

#include <memory>
#include <string>
#include <vector>

#include "../lib/tetgen.h"
#include "Femocs.h"

using namespace std;
namespace femocs {

/**
 * Class to create and handle FEM mesh in tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/
 */
class Mesh {
public:
    Mesh();

    void init_nodes(const int N);
    void init_faces(const int N);
    void init_elems(const int N);

    void init_facemarkers(const int N);
    void init_elemmarkers(const int N);
    void init_nodemarkers(const int N);

    void init_volumes(const int N);


    void add_node(const double x, const double y, const double z);
    void add_face(const int f1, const int f2, const int f3);
    void add_elem(const int e1, const int e2, const int e3, const int e4);

    void add_nodemarker(const int m);
    void add_facemarker(const int m);
    void add_elemmarker(const int m);

    void add_volume(const double V);

    void copy_statistics(shared_ptr<Mesh> mesh);

    void copy_nodes(shared_ptr<Mesh> mesh);
    void copy_faces(shared_ptr<Mesh> mesh, const int offset);
    void copy_elems(shared_ptr<Mesh> mesh, const int offset);

    void copy_nodemarkers(shared_ptr<Mesh> mesh);
    void copy_facemarkers(shared_ptr<Mesh> mesh);
    void copy_elemmarkers(shared_ptr<Mesh> mesh);

    const int getNnodes();
    const int getNelems();
    const int getNfaces();

    const int getNnodemarkers();
    const int getNfacemarkers();
    const int getNelemmarkers();

    const int getNvolumes();

    const int getNodemarker(const int i);
    const int getFacemarker(const int i);
    const int getElemmarker(const int i);
    double getVolume(const int i);

    int* getNodemarkers();
    vector<int>* getFacemarkers();
    vector<int>* getElemmarkers();

    double* getNodes();
    int* getFaces();
    int* getElems();

    const double getNode(int i, int xyz);
    const int getFace(int i, int node);
    const int getElem(int i, int node);

    void calc_volumes();
    void calc_volume_statistics();

    void recalc(const string cmd);
    void output(const string cmd);

    void write_faces(const string file_name);
    void write_elems(const string file_name);

    void transform_elemmarkers();

    tetgenio tetIO;
    /** Struct holding data about mesh statistics */
    struct Stat {
        double Vmin;
        double Vmax;
        double Vmedian;
        double Vaverage;
    };
    Stat stat;
private:
    tetgenbehavior tetgenbeh;

    int inodes;
    int ielems;
    int ifaces;

    int inodemarker;

    vector<int> nodemarkers;
    vector<int> facemarkers;
    vector<int> elemmarkers;
    vector<double> volumes;

    // Function to output mesh in .vtk format
//    void write_vtk(const string file_name, const int nnodes, const int ncells, const int nmarkers,
//            const REAL* nodes, const int* cells, const vector<int>* markers, const int celltype,
//            const int nnodes_per_cell);
    void write_vtk(const string file_name, const int nnodes, const int ncells, const int nnodemarkers, const int nmarkers,
            const REAL* nodes, const int* cells, const int* nodemarkers, const vector<int>* markers, const int celltype,
            const int nnodes_in_cell);
};

} /* namespace femocs */

#endif /* MESH_H_ */
