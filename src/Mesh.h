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
 * Class describing the centre of tetrahedral element
 */
class Centre {
public:
    Centre(double x, double y, double z);

    bool is_equal(Centre centre);
    double x;
    double y;
    double z;

private:
   // double r2;
};

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
    void init_centres(const int N);

    void set_facemarker(const int i, const int m);

    void add_node(const double x, const double y, const double z);
    void add_face(const int f1, const int f2, const int f3);
    void add_elem(const int e1, const int e2, const int e3, const int e4);

    void add_nodemarker(const int m);
    void add_facemarker(const int m);
    void add_elemmarker(const int m);

    void add_volume(const double V);
    void add_centre(const double x, const double y, const double z);

    void copy_statistics(shared_ptr<Mesh> mesh);

    void copy_nodes(shared_ptr<Mesh> mesh);
    void copy_faces(shared_ptr<Mesh> mesh, const int offset);
    void copy_elems(shared_ptr<Mesh> mesh, const int offset);

    void copy_nodemarkers(shared_ptr<Mesh> mesh);
    void copy_facemarkers(shared_ptr<Mesh> mesh);
    void copy_elemmarkers(shared_ptr<Mesh> mesh);

    const int get_n_nodes();
    const int get_n_elems();
    const int get_n_faces();

    const int get_n_nodemarkers();
    const int get_n_facemarkers();
    const int get_n_elemmarkers();

    const int get_n_volumes();

    const int get_nodemarker(const int i);
    const int get_facemarker(const int i);
    const int get_elemmarker(const int i);
    double get_volume(const int i);
    Centre get_centre(const int i);

    int* get_nodemarkers();
    vector<int>* get_facemarkers();
    vector<int>* get_elemmarkers();

    double* get_nodes();
    int* get_faces();
    int* get_elems();

    const double get_x(int i);
    const double get_y(int i);
    const double get_z(int i);
    const double get_node(int i, int xyz);
    const int get_face(int i, int node);
    const int get_elem(int i, int node);

    void calc_centres();
    void calc_volumes();
    void calc_volume_statistics();

    void recalc(const string cmd);
    void output(const string cmd);

    void write_nodes(const string file_name);
    void write_faces(const string file_name);
    void write_elems(const string file_name);

    bool centre_found_in(const int i, shared_ptr<Mesh> mesh);

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
    vector<Centre> centres;

    // Function to output mesh in .vtk format
//    void write_vtk(const string file_name, const int nnodes, const int ncells, const int nmarkers,
//            const REAL* nodes, const int* cells, const vector<int>* markers, const int celltype,
//            const int nnodes_per_cell);
    void write_vtk    (const string file_name, const int nnodes, const int ncells, const int nnodemarkers, const int nmarkers,
            const REAL* nodes, const int* cells, const int* nodemarkers, const vector<int>* markers, const int celltype,
            const int nnodes_in_cell);
};

} /* namespace femocs */

#endif /* MESH_H_ */
