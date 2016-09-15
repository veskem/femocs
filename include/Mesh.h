/*
 * Mesh.h
 *
 *  Created on: 1.2.2016
 *      Author: veske
 */

#ifndef MESH_H_
#define MESH_H_

#include "Macros.h"
#include "Primitives.h"
#include "AtomReader.h"
#include "Media.h"
#include "Medium.h"
#include "Tetgen.h"

using namespace std;
namespace femocs {

/**
 * Class to create and handle FEM mesh in tetgen format
 * http://wias-berlin.de/software/tetgen/1.5/
 */
class Mesh {
public:
    Mesh(const string& mesher);
    ~Mesh();

    /** Function to generate simple mesh that consists of one element */
    const void generate_simple();

    const void init_nodes(const int N);
    const void init_edges(const int N);
    const void init_faces(const int N);
    const void init_elems(const int N);

    const void init_nodemarkers(const int N);
    const void init_edgemarkers(const int N);
    const void init_facemarkers(const int N);
    const void init_elemmarkers(const int N);

    const void add_node(const double x, const double y, const double z);
    const void add_edge(const int n1, const int n2);
    const void add_face(const int n1, const int n2, const int n3);
    const void add_elem(const int n1, const int n2, const int n3, const int n4);
    const void add_node(const Point3 &point);
    const void add_edge(const SimpleEdge& edge);
    const void add_face(const SimpleFace& face);
    const void add_elem(const SimpleElement& elem);

    const void add_nodemarker(const int m);
    const void add_edgemarker(const int m);
    const void add_facemarker(const int m);
    const void add_elemmarker(const int m);

    const void copy_statistics(Mesh* mesh);
    const void copy_nodes(Mesh* mesh);
    const void copy_nodes(Mesh* mesh, const vector<bool> &mask);
    const void copy_edges(Mesh* mesh);
    const void copy_edges(const vector<SimpleEdge> &edges, const vector<bool> &mask);
    const void copy_faces(Mesh* mesh);
    const void copy_elems(Mesh* mesh);

    const void copy_nodemarkers(Mesh* mesh);
    const void copy_edgemarkers(Mesh* mesh);
    const void copy_facemarkers(Mesh* mesh);
    const void copy_elemmarkers(Mesh* mesh);

    const int get_n_nodes();
    const int get_n_edges();
    const int get_n_elems();
    const int get_n_faces();

    const int get_n_nodemarkers();
    const int get_n_edgemarkers();
    const int get_n_facemarkers();
    const int get_n_elemmarkers();

    const int get_n_areas();
    const int get_n_volumes();
    const int get_n_qualities();

    const int get_nodemarker(const int i);
    const int get_edgemarker(const int i);
    const int get_facemarker(const int i);
    const int get_elemmarker(const int i);

    const double get_area(const int i);
    const double get_volume(const int i);
    const double get_quality(const int i);

    const double* get_nodes();
    const int* get_edges();
    const int* get_faces();
    const int* get_elems();

    const vector<int>* get_nodemarkers();
    const vector<int>* get_edgemarkers();
    const vector<int>* get_facemarkers();
    const vector<int>* get_elemmarkers();

    const vector<int> get_elem_neighbours(const int i);

    const Vec3 get_vec(const int i);
    const Point3 get_node(const int i);
    const SimpleEdge get_simpleedge(const int i);
    const SimpleFace get_simpleface(const int i);
    const SimpleElement get_simpleelem(const int i);

    const Point3 get_face_centre(int i);
    const Point3 get_elem_centre(int i);

    const void set_nodemarker(const int node, const int m);
    const void set_edgemarker(const int edge, const int m);
    const void set_facemarker(const int face, const int m);
    const void set_elemmarker(const int elem, const int m);

    const void calc_qualities_byface();
    const void calc_qualities_byelem();
    const void calc_areas();
    const void calc_volumes();
    const void calc_statistics(const int i);
    const void calc_statistics();
    const void recalc();
    const void recalc(const string& cmd);
    const void recalc(const string& cmd1, const string& cmd2);
    const Medium to_medium();

    const void write_tetgen(const string file_name);
    const void write_nodes(const string file_name);
    const void write_edges(const string file_name);
    const void write_faces(const string file_name);
    const void write_elems(const string file_name);

    // Tetgen data structure
    tetgenio tetIOin;
    tetgenio tetIOout;

    /** Struct holding data about mesh statistics */
    struct Stat {
        double Vmin;
        double Vmax;
        double Vmedian;
        double Vaverage;
        int n_bulk;         //!< Number of nodes in bulk material
        int n_surface;      //!< Number of nodes on the surface of material
        int n_vacuum;       //!< Number of nodes in vacuum
        double xmin;    //!< Minimelem_in_vacuumum value of x-coordinate
        double xmax;    //!< Maximum value of x-coordinate
        double xmean;   //!< Average value of x-coordinate
        double ymin;    //!< Minimum value of y-coordinate
        double ymax;    //!< Maximum value of y-coordinate
        double ymean;   //!< Average value of y-coordinate
        double zmin;    //!< Minimum value of z-coordinate
        double zmax;    //!< Maximum value of z-coordinate
        double zmean;   //!< Average value of z-coordinate
    };

    /** Struct holding the indexes about nodes with known locations.
     * It's useful in finding the initially inserted nodes,
     * because when Tetgen adds nodes to the mesh, it adds them to the end of node list.
     */
    struct Indexes {
        int surf_start;
        int surf_end;
        int bulk_start;
        int bulk_end;
        int vacuum_start;
        int vacuum_end;
        int tetgen_start;
    };
    Stat stat;
    Indexes indxs;

    const int n_coordinates = 3;
    const int n_nodes_per_edge = 2;
    const int n_nodes_per_face = 3;
    const int n_nodes_per_elem = 4;
    const int n_edges_per_node = 3;
    const int n_edges_per_face = 3;
    const int n_edges_per_elem = 6;
    const int n_faces_per_elem = 4;
private:
    string mesher;
    int i_nodes;
    int i_edges;
    int i_faces;
    int i_elems;

    vector<int> nodemarkers;
    vector<int> edgemarkers;
    vector<int> facemarkers;
    vector<int> elemmarkers;

    vector<double> volumes;
    vector<double> areas;
    vector<double> qualities;

    const void init_statistics();

    double determinant(const double* n1, const double* n2, const double* n3, const double* n4, bool ones);
    const double determinant(const Vec3 &n1, const Vec3 &n2, const Vec3 &n3);
    const double calc_face_quality(const Vec3 &node1, const Vec3 &node2, const Vec3 &node3);
    const double calc_volume(const int i);
    const double calc_area(const int i);

    // Function to output mesh in .vtk format
    const void write_vtk(const string file_name, const int n_nodes, const int n_cells,
            const int n_markers, const REAL* nodes, const int* cells, const vector<int>* markers,
            const int celltype, const int n_nodes_in_cell);

    // Function to output nodes in .xyz format
    const void write_xyz(const string file_name);
};

} /* namespace femocs */

#endif /* MESH_H_ */
