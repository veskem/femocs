/*
 * TetgenCells.h
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#ifndef TETGENCELLS_H_
#define TETGENCELLS_H_

#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "Medium.h"

#include <fstream>

using namespace std;
namespace femocs {

/** Template class for holding finite element cells;
 * dim specifies the dimensionality of the cell - 1-node, 2-line, 3-triangle etc. */
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
    virtual void init(const int N) {
        require(N >= 0, "Invalid number of cells: " + to_string(N));
        i_cells = 0;
        *n_cells_w = N;
    }

    /** Append cell to the mesh */
    virtual void append(const SimpleCell<dim> &cell) {}

    /** Initialize marker appending */
    void init_markers(const int N) {
        require(N >= 0, "Invalid number of markers: " + to_string(N));
        markers.clear();
        markers.reserve(N);
    }

    /** Append cell marker */
    void append_marker(const int &m) {
        expect(get_n_markers() < markers.capacity(), "Allocated size of markers exceeded!");
        markers.push_back(m);
    }

    /** Set the number of available cells */
    void set_counter(const int n_cells) {
        require(n_cells >= 0, "Invalid number of cells: " + to_string(n_cells));
        i_cells = n_cells;
    }

    /** Copy the nodes from write buffer to read buffer */
    virtual void recalc() {
        *n_cells_r = *n_cells_w;
        i_cells = *n_cells_w;
    }

    /** Copy the cells from another mesh; mask can be used to copy only specified cells -
     * i-th true|false in mask means i-th cell is copied|skipped. */
    void copy(const TetgenCells<dim>& cells, const vector<bool>& mask={}) {
        const int n_cells = cells.size();

        // In case of empty or non-aligned mask, copy all the cells
        if (n_cells != mask.size()) {
            init(n_cells);
            for (int i = 0; i < n_cells; ++i)
                append(cells[i]);

        // In case of aligned mask, copy only the cells specified by the mask
        } else {
            const int n_mask = vector_sum(mask);
            init(n_mask);
            for (int i = 0; i < n_cells; ++i)
                if (mask[i])
                    append(cells[i]);
        }
    }

    /** Copy the cell markers from another mesh; mask can be used to copy only
     * specified markers - i-th true|false in mask means i-th marker is copied|skipped. */
    void copy_markers(const TetgenCells<dim>& cells, const vector<bool>& mask={}) {
        const int n_markers = cells.get_n_markers();

        // In case of empty or non-aligned mask, copy all the cell markers
        if (n_markers != mask.size()) {
            init_markers(n_markers);
            for (int i = 0; i < n_markers; ++i)
                append_marker(cells.get_marker(i));

        // In case of aligned mask, copy only the markers specified by the mask
        } else {
            const int n_mask = vector_sum(mask);
            init_markers(n_mask);
            for (int i = 0; i < n_markers; ++i)
                if (mask[i])
                    append_marker(cells.get_marker(i));
        }
    }

    /** Get number of cells in mesh */
    virtual int size() const { return *n_cells_r; }

    /** Get number of cell markers */
    int get_n_markers() const { return markers.size(); };

    /** Return i-th marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < get_n_markers(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    /** Get cell centroid coordinates */
    Point3 get_centroid(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));

        Point3 centroid(0.0);
        for (int v : get_cell(i))
            centroid += get_node(v);
        centroid *= 1.0 / dim;

        return centroid;
    }

    /** Return pointer to markers */
    vector<int>* get_markers() { return &markers; }

    /** Assign m-th marker */
    void set_marker(const int node, const int m) {
        require(node >= 0 && node < get_n_markers(), "Invalid index: " + to_string(node));
        markers[node] = m;
    }

    /** Function to write cell data to file */
    void write(const string &file_name) const {
        if (!MODES.WRITEFILE) return;

        string file_type = get_file_type(file_name);
        require(file_type == "vtk", "Unknown file type: " + file_type);

        int celltype = 0;
        if (dim == 1) celltype = 1;       // cell == vertex
        else if (dim == 2) celltype = 3;  // cell == line
        else if (dim == 3) celltype = 5;  // cell == triangle
        else if (dim == 4) celltype = 10; // cell == tetrahedron
        else if (dim == 8) celltype = 12; // cell == hexahedron

        write_vtk(file_name, celltype);
    }

    /** Accessor for accessing i-th cell */
    SimpleCell<dim> operator [](const size_t i) const { return get_cell(i); }

    /** Iterator to access the cells */
    typedef Iterator<TetgenCells, SimpleCell<dim>> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    const size_t DIM = dim;  //!< dimensionality of the cell; 1-node, 2-edge, 3-triangle etc

protected:
    const int n_coordinates = 3;  //!< number of spatial coordinates

    int* n_cells_r;      ///< number of readable cells in mesh data
    int* n_cells_w;      ///< number of writable cells in mesh data
    tetgenio* reads;     ///< mesh data that has been processed by Tetgen
    tetgenio* writes;    ///< mesh data that will be fed to Tetgen

    int i_cells;         ///< cell counter
    vector<int> markers; ///< cell markers

    /** Return i-th cell */
    virtual SimpleCell<dim> get_cell(const int i) const { return SimpleCell<dim>(); }

    /** Return number of readable nodes in the mesh */
    int get_n_nodes() const {
        return reads->numberofpoints;
    }

    /** Return i-th readable node from the mesh in point form */
    Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Point3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
    }
    
    /** Return i-th readable node from the mesh in vector form */
    Vec3 get_vec(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Vec3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
    }

    /** Output mesh in .vtk format */
    void write_vtk(const string &file_name, const int celltype) const {
        const int n_markers = get_n_markers();
        const int n_nodes = get_n_nodes();
        const int n_cells = size();

        expect(n_nodes > 0, "Zero nodes detected!");

        std::ofstream out(file_name.c_str());
        require(out, "Can't open a file " + file_name);

        out << fixed;

        out << "# vtk DataFile Version 3.0\n";
        out << "# TetgenCells data\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n\n";

        // Output the nodes
        out << "POINTS " << n_nodes << " double\n";
        for (size_t node = 0; node < n_nodes; ++node)
            out << get_node(node) << "\n";

        // Output the cells (tetrahedra, triangles, edges or vertices)
        out << "\nCELLS " << n_cells << " " << n_cells * (dim + 1) << "\n";
        for (size_t cl = 0; cl < n_cells; ++cl)
            out << dim << " " << get_cell(cl) << "\n";

        // Output cell types
        out << "\nCELL_TYPES " << n_cells << "\n";
        for (size_t cl = 0; cl < n_cells; ++cl)
            out << celltype << "\n";

        // Output cell markers
        if ((n_markers > 0) && (n_markers == n_cells)) {
            if (celltype == 1) out << "\nPOINT_DATA " << n_cells << "\n";
            else out << "\nCELL_DATA " << n_cells << "\n";
            out << "SCALARS Cell_markers int\nLOOKUP_TABLE default\n";
            for (size_t cl = 0; cl < n_cells; ++cl)
                out << get_marker(cl) << "\n";
        }
    }
};

/** Class for holding Tetgen nodes */
class TetgenNodes: public TetgenCells<1> {
public:
    TetgenNodes() : TetgenCells<1>() {}
    TetgenNodes(tetgenio *data) : TetgenCells<1>(data, &data->numberofpoints) {}
    TetgenNodes(tetgenio *read, tetgenio *write) :
        TetgenCells<1>(read, write, &read->numberofpoints, &write->numberofpoints) {}

    /** Modify the coordinates of i-th node */
    void set_node(const int i, const Point3 &point);

    /** Initialize node appending */
    void init(const int N);

    /** Append node to mesh */
    void append(const Point3 &point);

    /** Define accessor for accessing i-th node */
    Point3 operator [](const size_t i) const { return get_node(i); }

    /** Attach iterator */
    typedef Iterator<TetgenNodes, Point3> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Copy the nodes from write buffer to read buffer */
    void recalc();

    void copy(const TetgenNodes& nodes, const vector<bool>& mask={});

    /** Return the coordinates of i-th node as a 3D vector */
    Vec3 get_vec(const int i) const;

    /** Write node data to file */
    void write(const string &file_name) const;

    /** Save the locations of the initially added nodes */
    void save_indices(const int n_surf, const int n_bulk, const int n_vacuum);

    /** Store the locations of different kinds of nodes that were produced while splitting tetrahedra into hexahedra */
    void save_hex_indices(const vector<int>& n_nodes);

    /** Calculate statistics about nodes */
    void calc_statistics();

    vector<dealii::Point<3>> export_dealii() const;

    /** Struct holding the indexes about nodes with known locations.
     * It's useful for finding the initially inserted nodes,
     * because when Tetgen adds nodes to the mesh, it adds them to the end of the pointlist. */
    struct Indexes {
        int surf_start, surf_end;
        int bulk_start, bulk_end;
        int vacuum_start, vacuum_end;
        int tetgen_start, tetgen_end;

        int tetnode_start, tetnode_end;
        int midface_start, midface_end;
        int midedge_start, midedge_end;
        int midtet_start, midtet_end;

        int size() const { return (&midtet_end)-(&surf_start) + 1; };  ///< number of values
    } indxs;

    /** Struct holding statistics about nodes */
    struct Stat {
        int n_tetnode;
        int n_midface;
        int n_midedge;
        int n_midtet;

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
    SimpleCell<1> get_cell(const int i) const;

    void copy_statistics(const TetgenNodes& nodes);

    /** Write node data to .xyz file */
    void write_xyz(const string &file_name) const;

    /** Initialize statistics about nodes */
    void init_statistics();
};

/** Class for holding Tetgen line edges */
class TetgenEdges: public TetgenCells<2> {
public:
    TetgenEdges() : TetgenCells<2>() {}
    TetgenEdges(tetgenio *data) : TetgenCells<2> (data, &data->numberofedges) {}
    TetgenEdges(tetgenio *read, tetgenio *write)
        : TetgenCells<2>(read, write, &read->numberofedges, &write->numberofedges) {}

    /** Initialize edge appending */
    void init(const int N);

    /** Append edge to mesh */
    void append(const SimpleCell<2> &cell);

    /** Copy the nodes from write buffer to read buffer */
    void recalc();

private:
    /** Return i-th edge */
    SimpleCell<2> get_cell(const int i) const;
};

/** Class for holding Tetgen triangular faces */
class TetgenFaces: public TetgenCells<3> {
public:
    TetgenFaces() : TetgenCells<3>() {}
    TetgenFaces(tetgenio *data) : TetgenCells<3> (data, &data->numberoftrifaces) {}
    TetgenFaces(tetgenio *read, tetgenio *write)
        : TetgenCells<3>(read, write, &read->numberoftrifaces, &write->numberoftrifaces) {}

    /** Initialize face appending */
    void init(const int N);

    /** Append face to mesh */
    void append(const SimpleCell<3> &cell);

    /** Copy the nodes from write buffer to read buffer */
    void recalc();
    
    /** Delete the faces on the sides of simulation cell */
    void clean_sides(const Medium::Sizes& stat);

    /** Return the normal of i-th triangle */
    Vec3 get_norm(const int i) const;
    
    /** Return the area of i-th triangle */
    double get_area(const int i) const;

private:
    vector<double> areas;   ///< areas of triangles
    vector<Vec3> norms;     ///< norms of triangles

    /** Return i-th face */
    SimpleCell<3> get_cell(const int i) const;

    /** Calculate the norms and areas for all the triangles */
    void calc_norms();
};

/** Class for holding Tetgen tetrahedral elements */
class TetgenElements: public TetgenCells<4> {
public:
    TetgenElements() : TetgenCells<4>() {}
    TetgenElements(tetgenio *data) : TetgenCells<4> (data, &data->numberoftetrahedra) {}
    TetgenElements(tetgenio *read, tetgenio *write) :
        TetgenCells<4>(read, write, &read->numberoftetrahedra, &write->numberoftetrahedra) {}

    /** Get indices of neighbouring elements of i-th element */
    vector<int> get_neighbours(const int i) const;

    /** Initialize element appending */
    void init(const int N);

    /** Append element to mesh */
    void append(const SimpleCell<4> &cell);

    /** Copy the nodes from write buffer to read buffer */
    void recalc();

private:
    /** Return i-th element */
    SimpleCell<4> get_cell(const int i) const;
};

class Hexahedra: public TetgenCells<8> {
public:
    /** SimpleCells constructors */
    Hexahedra() : TetgenCells<8>() {}
    Hexahedra(tetgenio *data) : TetgenCells<8> (data, &data->numberofvcells) {}
    Hexahedra(tetgenio *read, tetgenio *write) :
        TetgenCells<8>(read, write, &read->numberofvcells, &write->numberofvcells) {}

    /** Initialize cell appending */
    void init(const int N);

    /** Append cell to the mesh */
    void append(const SimpleCell<8> &cell);

    /** Get number of cells in mesh */
    int size() const;

    vector<dealii::CellData<3>> export_dealii() const;

protected:
    vector<SimpleHex> hexs;

    /** Return i-th hexahedron */
    SimpleCell<8> get_cell(const int i) const;
};

/** Virtual class for holding data that is common to Voronoi cell and Voronoi face */
class Voronoi {
public:
    Voronoi() : data(NULL), id(-1) {}
    Voronoi(tetgenio* data, const int i) : data(data), id(i) {}
    virtual ~Voronoi() {}

    /** Get number of faces that make the cell */
    virtual int size() const { return 0; }

    int id;

protected:
    tetgenio* data;     ///< mesh data that has been processed by Tetgen
};

/** Class for accessing the Voronoi face data */
class VoronoiFace : public Voronoi {
public:
    VoronoiFace() : Voronoi() {}
    VoronoiFace(tetgenio* data, const int i, const int calc=0) : Voronoi(data, i) {
        if (calc) calc_nodes();
    }
    ~VoronoiFace() {}

    /** Get number of nodes that make the face */
    int size() const { return data->vfacetlist[id].elist[0]; }

    /** Get the norm vector of the face */
    Vec3 norm();

    /** Get the area of the face.
     * See the theory in http://geomalgorithms.com/a01-_area.html#3D%20Polygons */
    double area();

    /** Accessor for accessing the i-th node */
    int operator [](size_t i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return nodes[i];
    }

    /** Iterator to access the nodes of faces */
    typedef Iterator<VoronoiFace, int> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Stream for printing the nodes to file or console */
    friend std::ostream& operator <<(std::ostream &os, const VoronoiFace &vf) {
        for (int n : vf.nodes) os << n << ' ';
        return os;
    }

    /** Transform the node data from tetgenio into easily accessible form */
    void calc_nodes();

private:
    vector<Vec3> verts;  ///< coordinates of the face vertices
    vector<int> nodes;   ///< indices of the face nodes
};

/** Class for accessing the Voronoi cell data */
class VoronoiCell : public Voronoi {
public:
    VoronoiCell() : Voronoi() {}
    VoronoiCell(tetgenio* data, const int i, const int calc=0) : Voronoi(data, i) {}
    ~VoronoiCell() {}

    /** Get number of faces that make the cell */
    int size() const { return data->vcelllist[id][0]; }

    /** Accessor for accessing the i-th face */
    VoronoiFace operator [](size_t i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return VoronoiFace(data, data->vcelllist[id][i+1]);
    }

    /** Iterator to access the cell faces */
    typedef Iterator<VoronoiCell, VoronoiFace> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

private:
};

/** Virtual class for holding data that is common to Voronoi cells and Voronoi faces */
template<typename T>
class Voronois {
public:
    /** Voronoi constructors */
    Voronois() : tetio(NULL), _n_cells(NULL) {}

    Voronois(tetgenio* data, int* n_cells) : tetio(data), _n_cells(n_cells) {}

    /** Voronoi destructor */
    virtual ~Voronois() {};

    /** Initialize marker appending */
    void init_markers(const int N) {
        require(N >= 0, "Invalid number of markers: " + to_string(N));
        markers.clear();
        markers.reserve(N);
    }

    /** Append cell marker */
    void append_marker(const int &m) {
        expect(get_n_markers() < markers.capacity(), "Allocated size of markers exceeded!");
        markers.push_back(m);
    }

    /** Get number of cells in mesh */
    int size() const { return *_n_cells; }

    /** Get number of cell markers */
    int get_n_markers() const { return markers.size(); };

    /** Return i-th marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < get_n_markers(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    /** Return pointer to markers */
    vector<int>* get_markers() { return &markers; }

    /** Assign m-th marker */
    void set_marker(const int node, const int m) {
        require(node >= 0 && node < get_n_markers(), "Invalid index: " + to_string(node));
        markers[node] = m;
    }

    /** Function to write cell data to file */
    void write(const string &file_name) const {
        if (!MODES.WRITEFILE) return;

        string file_type = get_file_type(file_name);
        if (file_type == "vtk")
            write_vtk(file_name);
        else
            require(false, "Unknown file type: " + file_type);
    }

    /** Accessor for accessing the i-th cell */
    T operator [](size_t i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return T(tetio, i);
    }

    /** Iterator to access the cells */
    typedef Iterator<Voronois, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

protected:
    const int n_coordinates = 3;  ///< number of spatial coordinates
    const int celltype = 7;       ///< vtk cell = polygon

    int* _n_cells;        ///< number of cells in mesh data
    tetgenio* tetio;     ///< mesh data that has been processed by Tetgen
    vector<int> markers; ///< cell markers

    /** Return number of readable nodes in the mesh */
    int get_n_nodes() const { return tetio->numberofvpoints; }

    /** Return i-th node from the voronoi mesh */
    Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Point3(tetio->vpointlist[n+0], tetio->vpointlist[n+1], tetio->vpointlist[n+2]);
    }

    Vec3 get_vec(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + to_string(i));
        const int n = n_coordinates * i;
        return Vec3(tetio->vpointlist[n+0], tetio->vpointlist[n+1], tetio->vpointlist[n+2]);
    }

    void write_vtk(const string &file_name) const {
        const int n_markers = get_n_markers();
        const int n_nodes = get_n_nodes();

        expect(n_nodes > 0, "Zero nodes detected!");

        std::ofstream out(file_name.c_str());
        require(out, "Can't open a file " + file_name);

        out << fixed;

        out << "# vtk DataFile Version 3.0\n";
        out << "Voronoi cells\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n\n";

        // Output the nodes
        out << "POINTS " << n_nodes << " double\n";
        for (size_t node = 0; node < n_nodes; ++node)
            out << get_node(node) << "\n";

        // Output the cells and associated data
        write_cells(out);
    }

    virtual void write_cells(ofstream& out) const {}
};

/** Class for accessing the faces of Voronoi cells */
class VoronoiFaces : public Voronois<VoronoiFace> {
public:
    /** VoronoiCells constructors */
    VoronoiFaces() : Voronois<VoronoiFace>() {}
    VoronoiFaces(tetgenio* data) : Voronois<VoronoiFace>(data, &data->numberofvfacets) {}

    void write_cells(ofstream& out) const;
};

/** Class for accessing the Voronoi cells */
class VoronoiCells : public Voronois<VoronoiCell> {
public:
    /** VoronoiCells constructors */
    VoronoiCells() : Voronois<VoronoiCell>() {}
    VoronoiCells(tetgenio* data) : Voronois<VoronoiCell>(data, &data->numberofvcells) {}

    void write_cells(ofstream& out) const;
};

} /* namespace femocs */

#endif /* TETGENCELLS_H_ */
