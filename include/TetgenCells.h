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
#include <fstream>

using namespace std;
namespace femocs {

/** Template class for holding finite element cells */
template<size_t dim>
class TetgenCells {
public:
    /** SimpleCells constructors */
    TetgenCells() : i_cells(0), reads(NULL), writes(NULL), n_cells_r(NULL), n_cells_w(NULL) {}

    /** Constructor for the case when read and write data go into same place */
    TetgenCells(tetgenio* data, int* n_cells) :
        i_cells(0), reads(data), writes(data), n_cells_r(n_cells), n_cells_w(n_cells) {}

    /** Constructor for the case when read and write data go into different places */
    TetgenCells(tetgenio* reads, tetgenio* writes, int* n_cells_r, int* n_cells_w) :
        i_cells(0), reads(reads), writes(writes), n_cells_r(n_cells_r), n_cells_w(n_cells_w) {}

    /** SimpleCells destructor */
    virtual ~TetgenCells() {};

    /** Initialize cell appending */
    const void init(const int N) {
        require(N > 0, "Invalid number of cells: " + to_string(N));
        i_cells = 0;
        *n_cells_w = N;
    }

    /** Append cell to the mesh */
    virtual const void append(const SimpleCell<dim> &cell) {}

    /** Initialize marker appending */
    const void init_markers(const int N) {
        require(N > 0, "Invalid number of markers: " + to_string(N));
        markers.reserve(N);
    }

    /** Append cell marker */
    const void append_marker(const int &m) {
        expect(get_n_markers() < markers.capacity(), "Allocated size of markers exceeded!");
        markers.push_back(m);
    }

    /** Set the number of available cells */
    void set_counter(const int n_cells) {
        require(n_cells >= 0, "Invalid number of cells: " + to_string(n_cells));
        i_cells = n_cells;
    }

    /** Copy the cells (nodes, edges, faces or elements) from another mesh;
     * mask can be used to copy only the specified cells -
     * i-th true in means i-th cell is copied and false means it's skipped. */
    const void copy(const TetgenCells<dim>& cells, const vector<bool>& mask={}) {
        const int n_nodes = cells.get_n_cells();

        // In case of empty or non-aligned mask, copy all the cells
        if (n_nodes != mask.size()) {
            init(n_nodes);
            for (int i = 0; i < n_nodes; ++i)
                append(cells[i]);
            i_cells = n_nodes;

        // In case of aligned mask, copy only the cells specified by the mask
        } else {
            const int n_mask = vector_sum(mask);
            init(n_mask);
            for (int i = 0; i < n_nodes; ++i)
                if (mask[i])
                    append(cells[i]);
            i_cells = n_mask;
        }
    }

    const void copy_markers(const TetgenCells<dim>& cells, const vector<bool>& mask={}) {
        const int n_markers = cells.get_n_markers();

        // In case of empty or non-aligned mask, copy all the cell markers
        if (n_markers != mask.size()) {
            init_markers(n_markers);
            for (int i = 0; i < n_markers; ++i)
                append_marker(cells->get_marker(i));

        // In case of aligned mask, copy only the markers specified by the mask
        } else {
            const int n_mask = vector_sum(mask);
            init_markers(n_mask);
            for (int i = 0; i < n_markers; ++i)
                if (mask[i])
                    append_marker(cells->get_marker(i));
        }
    }

    /** Return number of nodes in mesh */
    const int get_n_nodes() const { return reads->numberofpoints; }

    /** Get number of cells in mesh */
    const int get_n_cells() const { return *n_cells_r; }

    /** Get number of cell markers */
    const int get_n_markers() const { return markers.size(); };

    /** Return i-th node from mesh */
    const Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Point3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
    }

    /** Return i-th marker */
    const int get_marker(const int i) const {
        require(i >= 0 && i < get_n_markers(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    /** Get cell centroid coordinates */
    const Point3 get_centroid(const int i) const {
        require(i >= 0 && i < get_n_cells(), "Invalid index: " + to_string(i));

        Point3 centroid(0.0);
        for (int v : get_cell(i))
            centroid += get_node(v);
        centroid *= 1.0 / dim;

        return centroid;
    }

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

    /** Define accessor for accessing i-th cell */
    const SimpleCell<dim> operator [](size_t i) const { return get_cell(i); }

    /** Attach iterator */
    typedef Iterator<TetgenCells, SimpleCell<dim>> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, get_n_cells()); }

    const size_t DIM = dim;

    /** Pointers to Tetgen data */
    const int* n_cells_r;
    int* n_cells_w;
    tetgenio* reads;
    tetgenio* writes;

protected:
    const int n_coordinates = 3;



    /** Internal variables */
    int i_cells;
    vector<int> markers;

    /** Return i-th cell */
    virtual const SimpleCell<dim> get_cell(const int i) const { return SimpleCell<dim>(); }

    /** Output mesh in .vtk format */
    const void write_vtk(const string &file_name, const int celltype) {
        const int n_markers = get_n_markers();
        const int n_nodes = get_n_nodes();
        const int n_cells = get_n_cells();

        expect(n_nodes > 0, "Zero nodes detected!");

        std::ofstream out(file_name.c_str());
        require(out, "Can't open a file " + file_name);

        out.setf(std::ios::scientific);
        out.precision(8);

        out << "# vtk DataFile Version 3.0\n";
        out << "# TetgenCells data\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n\n";

        // Output the nodes
        out << "POINTS " << n_nodes << " double\n";
        for (size_t ver = 0; ver < n_nodes; ++ver)
            out << get_node(ver) << "\n";

        // Output the cells (tetrahedra, triangles, edges or vertices)
        out << "\nCELLS " << n_cells << " " << n_cells * (dim + 1) << "\n";
        for (size_t el = 0; el < n_cells; ++el)
            out << dim << " " << get_cell(el) << "\n";

        // Output cell types
        out << "\nCELL_TYPES " << n_cells << "\n";
        for (size_t el = 0; el < n_cells; ++el)
            out << celltype << "\n";

        // Output cell markers
        if ((n_markers > 0) && (n_markers == n_cells)) {
            out << "\nCELL_DATA " << n_cells << "\n";
            out << "SCALARS Cell_markers int\nLOOKUP_TABLE default\n";
            for (size_t el = 0; el < n_cells; ++el)
                out << get_marker(el) << "\n";
        }
    }
};

/** Class for holding simple nodes */
class TetgenNodes: public TetgenCells<1> {
public:
    TetgenNodes() : TetgenCells<1>() {}
    TetgenNodes(tetgenio *data) : TetgenCells<1>(data, &data->numberofpoints) {}
    TetgenNodes(tetgenio *read, tetgenio *write) :
        TetgenCells<1>(read, write, &read->numberofpoints, &write->numberofpoints) {}

    /** Initialize node appending */
    const void init(const int N) {
        TetgenCells::init(N);
        writes->pointlist = new double[n_coordinates * N];
    }

    /** Append node to mesh */
    const void append(const Point3 &point) {
        require(i_cells < *n_cells_w, "Allocated size of nodes exceeded!");
        int i = n_coordinates * i_cells;
        for (double node : point)
            writes->pointlist[i++] = node;
        i_cells++;
    }

    /** Return the coordinates of i-th node as a 3D vector */
    const Vec3 get_vec(const int i) const;

    /** Write node data to file */
    const void write(const string &file_name);

    const void save_indices(const int n_surf, const int n_bulk, const int n_vacuum);

    const void init_statistics();

    const void calc_statistics();

    /** Struct holding the indexes about nodes with known locations.
     * It's useful in finding the initially inserted nodes,
     * because when Tetgen adds nodes to the mesh, it adds them to the end of node list.
     */
    struct Indexes {
        int surf_start, surf_end;
        int bulk_start, bulk_end;
        int vacuum_start, vacuum_end;
        int tetgen_start;
    } indxs;

    /** Struct holding statistics about nodes */
    struct Stat {
        int n_bulk;     //!< Number of nodes in bulk material
        int n_surface;  //!< Number of nodes on the surface of material
        int n_vacuum;   //!< Number of nodes in vacuum

        double xmin;    //!< Minimum value of x-coordinate
        double xmax;    //!< Maximum value of x-coordinate
        double ymin;    //!< Minimum value of y-coordinate
        double ymax;    //!< Maximum value of y-coordinate
        double zmin;    //!< Minimum value of z-coordinate
        double zmax;    //!< Maximum value of z-coordinate
    } stat;

private:
    /** Return index of i-th node */
    const SimpleCell<1> get_cell(const int i) const;

    /** Write node data to .xyz file */
    const void write_xyz(const string &file_name);
};

/** Class for holding Tetgen line edges */
class TetgenEdges: public TetgenCells<2> {
public:
    TetgenEdges() : TetgenCells<2>() {}
    TetgenEdges(tetgenio *data) : TetgenCells<2> (data, &data->numberofedges) {}
    TetgenEdges(tetgenio *read, tetgenio *write)
        : TetgenCells<2>(read, write, &read->numberofedges, &write->numberofedges) {}

    /** Initialize edge appending */
    const void init(const int N) {
        TetgenCells::init(N);
        writes->edgelist = new int[DIM * N];
    }

    /** Append edge to mesh */
    const void append(const SimpleEdge &cell) {
        require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
        int i = DIM * i_cells;
        for (unsigned int node : cell)
            writes->edgelist[i++] = node;
        i_cells++;
    }

private:
    /** Return i-th edge */
    const SimpleCell<2> get_cell(const int i) const;
};

/** Class for holding Tetgen triangular faces */
class TetgenFaces: public TetgenCells<3> {
public:
    TetgenFaces() : TetgenCells<3>() {}
    TetgenFaces(tetgenio *data) : TetgenCells<3> (data, &data->numberoftrifaces) {}
    TetgenFaces(tetgenio *read, tetgenio *write)
        : TetgenCells<3>(read, write, &read->numberoftrifaces, &write->numberoftrifaces) {}

    /** Initialize face appending */
    const void init(const int N) {
        TetgenCells::init(N);
        writes->trifacelist = new int[DIM * N];
    }

    /** Append face to mesh */
    const void append(const SimpleFace &cell) {
        require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
        int i = DIM * i_cells;
        for (unsigned int node : cell)
            writes->trifacelist[i++] = node;
        i_cells++;
    }

private:
    /** Return i-th face */
    const SimpleCell<3> get_cell(const int i) const;
};

/** Class for holding Tetgen tetrahedral elements */
class TetgenElements: public TetgenCells<4> {
public:
    TetgenElements() : TetgenCells<4>() {}
    TetgenElements(tetgenio *data) : TetgenCells<4> (data, &data->numberoftetrahedra) {}
    TetgenElements(tetgenio *read, tetgenio *write) :
        TetgenCells<4>(read, write, &read->numberoftetrahedra, &write->numberoftetrahedra) {}

    /** Get indices of neighbouring elements of i-th element */
    const vector<int> get_neighbours(const int i) const;

    /** Initialize element appending */
    const void init(const int N) {
        TetgenCells::init(N);
        writes->tetrahedronlist = new int[DIM * N];
    }

    /** Append element to mesh */
    const void append(const SimpleElement &cell) {
        require(i_cells < *n_cells_w, "Allocated size of cells exceeded!");
        int i = DIM * i_cells;
        for (unsigned int node : cell)
            writes->tetrahedronlist[i++] = node;
        i_cells++;
    }
private:
    /** Return i-th element */
    const SimpleCell<4> get_cell(const int i) const;
};

} /* namespace femocs */

#endif /* SIMPLECELL_H_ */
