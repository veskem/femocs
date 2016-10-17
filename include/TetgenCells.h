/*
 * SimpleCell.h
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#ifndef SIMPLECELL_H_
#define SIMPLECELL_H_

#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"

using namespace std;
namespace femocs {

/** Template class for holding finite element cells */
template<size_t dim>
class TetgenCells {
public:
    /** SimpleCells constructors */
    TetgenCells() :
        i_cells(0),
        points_r(NULL), n_points_r(NULL),
        nodes_r(NULL), nodes_w(NULL),
        n_cells_r(NULL), n_cells_w(NULL) {}

    /** Constructor for the case when read and write data go into same place */
    TetgenCells(double* points, int* n_points, int* nodes, int* n_cells) :
        i_cells(0),
        points_r(points), n_points_r(n_points),
        nodes_r(nodes), nodes_w(nodes),
        n_cells_r(n_cells), n_cells_w(n_cells) {}

    /** Constructor for the case when read and write data go into different places */
    TetgenCells(double* points_r, int* n_points_r, int* nodes_r, int* nodes_w, int* n_cells_r, int* n_cells_w) :
        i_cells(0),
        points_r(points_r), n_points_r(n_points_r),
        nodes_r(nodes_r), nodes_w(nodes_w),
        n_cells_r(n_cells_r), n_cells_w(n_cells_w) {}

    /** SimpleCells destructor */
    virtual ~TetgenCells() {};

    /** Initialize cell appending */
    const void init_cells(const int N) {
        require(N > 0, "Invalid number of cells: " + to_string(N));
        i_cells = 0;
        n_cells_r = N;
        nodes_r = new int[dim * N];
    }

    /** Initialize marker appending */
    const void init_markers(const int N) {
        require(N > 0, "Invalid number of markers: " + to_string(N));
        markers.reserve(N);
    }

    /** Append cell to mesh */
    const void add_cell(const SimpleCell<dim> &cell) {
        require(i_cells < *n_cells_r, "Allocated size of cells exceeded!");
        require(vector_sum(cell >= 0) == dim, "Invalid cell: " + cell.to_str());

        int i = dim * i_cells;
        for (unsigned int node : cell)
            nodes_r[i++] = node;
        i_cells++;
    }

    /** Append cell marker */
    const void add_marker(const int &m) {
        expect(get_n_markers() < markers.capacity(), "Allocated size of markers exceeded!");
        markers.push_back(m);
    }

//    const void copy_cells(Mesh* mesh, const vector<bool> &mask={});
//    const void copy_markers(Mesh* mesh, const vector<bool> &mask={});

    /** Return number of nodes in mesh */
    const int get_n_nodes() const { return *n_points_r; }

    /** Get number of cells in mesh */
    const int get_n_cells() const { return *n_cells_r; }

    /** Get number of cell markers */
    const int get_n_markers() const { return markers.size(); };

    /** Return i-th node from mesh */
    const Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Point3(points_r[n+0], points_r[n+1], points_r[n+2]);
    }

    /** Return i-th cell */
    virtual const SimpleCell<dim> get_cell(const int i) { return SimpleCell<dim>(); }

    /** Return i-th marker */
    const int get_marker(const int i) const {
        require(i >= 0 && i < get_n_markers(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    /** Get cell centroid coordinates */
    const Point3 get_centroid(const int i) const {
        require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));

        Point3 centroid;
        for (int v : get_cell(i))
            centroid += get_node(v);
        centroid *= 1.0 / dim;

        return centroid;
    }

    /** Return pointer to cells */
    const int* get_cells() { return nodes_r; }

    /** Return pointer to markers */
    const vector<int>* get_markers() { return &markers; }

    /** Assign m-th marker */
    const void set_marker(const int node, const int m) {
        require(node >= 0 && node < get_n_markers(), "Invalid index: " + to_string(node));
        markers[node] = m;
    }

    /** Function to write cell data to file */
    const void write(const string &file_name) {
    #if not DEBUGMODE
        return;
    #endif
        string file_type = get_file_type(file_name);
        require(file_type == "vtk", "Unknown file type: " + file_type);

        int celltype = 0;
        if (dim == 1) celltype = 1;       // cell == vertex
        else if (dim == 2) celltype = 3;  // cell == line
        else if (dim == 3) celltype = 5;  // cell == triangle
        else if (dim == 4) celltype = 10; // cell == tetrahedron

        write_vtk(file_name, celltype);
    }

protected:
    const int n_coordinates = 3;
    const size_t DIM = dim;

    /** Internal variables */
    int i_cells;
    vector<int> markers;

    /** Pointers to Tetgen data */
    const double* points_r;
    const int* n_points_r;
    const int* nodes_r;
    int* nodes_w;
    const int* n_cells_r;
    int* n_cells_w;

    /** Output mesh in .vtk format */
    const void write_vtk(const string &file_name, const int celltype) {
        const int n_markers = get_n_markers();
        const int n_nodes = get_n_nodes();
        const int n_cells = *n_cells_r;

        char file_name_char[1024];
        strcpy(file_name_char, file_name.c_str());

        expect(n_nodes > 0, "Zero nodes detected!");

        FILE *out_file;
        out_file = fopen(file_name_char, "w");
        require(out_file != (FILE*) NULL, "Can't open a file " + file_name);

        fprintf(out_file, "# vtk DataFile Version 3.0\n");
        fprintf(out_file, "# Unstructured grid\n");
        fprintf(out_file, "ASCII\n"); // another option is BINARY
        fprintf(out_file, "DATASET UNSTRUCTURED_GRID\n\n");

        // Output the nodes
        if (n_nodes > 0) {
            fprintf(out_file, "POINTS %d double\n", n_nodes);
            for (int i = 0; i < n_nodes; ++i)
                fprintf(out_file, "%s\n", get_node(i).to_str().c_str());
            fprintf(out_file, "\n");
        }

        // Output the cells (tetrahedra, triangles, edges or vertices)
        if (n_cells > 0) {
            fprintf(out_file, "CELLS %d %d\n", n_cells, n_cells * (dim + 1));
            for (int i = 0; i < n_cells; ++i)
                fprintf( out_file, "%d %s\n", (dim, get_cell(i).to_str().c_str()) );
            fprintf(out_file, "\n");
        }

        // Output the types of cells
        if (n_cells_r > 0) {
            fprintf(out_file, "CELL_TYPES %d\n", n_cells);
            for (int i = 0; i < n_cells; ++i)
                fprintf(out_file, "%d \n", celltype);
            fprintf(out_file, "\n\n");
        }

        // Output cell markers
        if ((n_markers > 0) && (n_markers == n_cells)) {
            fprintf(out_file, "CELL_DATA %d\n", n_markers);
            fprintf(out_file, "SCALARS Cell_markers int\n");
            fprintf(out_file, "LOOKUP_TABLE default\n");
            for (int i = 0; i < n_markers; ++i)
                fprintf(out_file, "%d\n", get_marker(i));
            fprintf(out_file, "\n");
        }

        fclose(out_file);
    }
};

/** Class for holding simple nodes */
class TetgenNodes: public TetgenCells<1> {
public:
    TetgenNodes() : TetgenCells<1>(), points_w(NULL), n_points_w(NULL) {}
    TetgenNodes(tetgenio *data) :
        TetgenCells<1>(data->pointlist, &data->numberofpoints, NULL, &data->numberofpoints),
        points_w(data->pointlist), n_points_w(&data->numberofpoints) {}

    TetgenNodes(tetgenio *read, tetgenio *write) :
        TetgenCells<1>(read->pointlist, &read->numberofpoints, NULL, NULL, &read->numberofpoints, NULL),
        points_w(write->pointlist), n_points_w(&write->numberofpoints) {}

    /** Return index of i-th node */
    const SimpleNode get_cell(const int i) const;

    /** init_cells does not have meaning for nodes, so it's disabled */
    const void init_cells(const int N) {
        require(false, "Invalid usage of SimpleNodes!");
    }

    /** Initialize node appending */
    const void init_nodes(const int N);

    /** add_cell does not have meaning for nodes, so it's disabled */
    const void add_cell(const SimpleCell<1> &cell) {
        require(false, "Invalid usage of SimpleNodes!");
    }

    /** Append node to mesh */
    const void add_node(const Point3 &point);

    /** Return the coordinates of i-th node as a 3D vector */
    const Vec3 get_vec(const int i) const;

    /** Write node data to file */
    const void write(const string &file_name);

private:
    /** Pointers to Tetgen data */
    double* points_w;
    int* n_points_w;

    /** Write node data to .xyz file */
    const void write_xyz(const string &file_name);
};

/** Class for holding Tetgen line edges */
class TetgenEdges: public TetgenCells<2> {
public:
    TetgenEdges() : TetgenCells<2>() {}
    TetgenEdges(tetgenio *data) : TetgenCells<2> (data->pointlist, &data->numberofpoints, data->edgelist, &data->numberofedges) {}
    TetgenEdges(tetgenio *read, tetgenio *write)
        : TetgenCells<2>(read->pointlist, &read->numberofpoints, read->edgelist, write->edgelist, &read->numberofedges, &write->numberofedges) {}

    /** Return i-th edge */
    const SimpleEdge get_cell(const int i) const;
};

/** Class for holding Tetgen triangular faces */
class TetgenFaces: public TetgenCells<3> {
public:
    TetgenFaces() : TetgenCells<3>() {}
    TetgenFaces(tetgenio *data) : TetgenCells<3> (data->pointlist, &data->numberofpoints, data->edgelist, &data->numberofedges) {}
    TetgenFaces(tetgenio *read, tetgenio *write)
        : TetgenCells<3>(read->pointlist, &read->numberofpoints, read->edgelist, write->edgelist, &read->numberofedges, &write->numberofedges) {}

    /** Return i-th face */
    const SimpleFace get_cell(const int i) const;
};

/** Class for holding Tetgen tetrahedral elements */
class TetgenElements: public TetgenCells<4> {
public:
    TetgenElements() : TetgenCells<4>(), neighborlist(NULL) {}

    TetgenElements(tetgenio *data) :
        TetgenCells<4> (data->pointlist, &data->numberofpoints, data->edgelist, &data->numberofedges),
        neighborlist(data->neighborlist) {}

    TetgenElements(tetgenio *read, tetgenio *write) :
        TetgenCells<4>(read->pointlist, &read->numberofpoints, read->edgelist, write->edgelist, &read->numberofedges, &write->numberofedges),
        neighborlist(read->neighborlist) {}

    /** Return i-th element */
    const SimpleElement get_cell(const int i) const;

    /** Get indices of neighbouring elements of i-th element */
    const vector<int> get_neighbours(const int i) const;

private:
    const int* neighborlist;
};

} /* namespace femocs */

#endif /* SIMPLECELL_H_ */