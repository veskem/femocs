/*
 * TetgenCells.h
 *
 *  Created on: 27.9.2016
 *      Author: veske
 */

#ifndef TETGENCELLS_H_
#define TETGENCELLS_H_

#include "FileWriter.h"
#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "Medium.h"

#include <fstream>

using namespace std;
namespace femocs {

/** Template class for holding finite element cells;
 * dim specifies the dimensionality of the cell - 1-node, 2-line, 3-triangle etc. */
template<int dim>
class TetgenCells: public FileWriter {
public:

    /** Default constructor that creates empty instance */
    TetgenCells() :
    n_cells_r(NULL), n_cells_w(NULL), reads(NULL), writes(NULL), i_cells(0) {}

    /** Constructor for the case when read and write data go into same place */
    TetgenCells(tetgenio* data, int* n_cells) :
        n_cells_r(n_cells), n_cells_w(n_cells), reads(data), writes(data), i_cells(0) {}

    /** Constructor for the case when read and write data go into different places */
    TetgenCells(tetgenio* reads, tetgenio* writes, int* n_cells_r, int* n_cells_w) :
        n_cells_r(n_cells_r), n_cells_w(n_cells_w), reads(reads), writes(writes), i_cells(0) {}

    /** SimpleCells destructor */
    virtual ~TetgenCells() {};

    /** Initialize cell appending */
    virtual void init(const int N) {
        require(N >= 0, "Invalid number of cells: " + d2s(N));
        i_cells = 0;
        *n_cells_w = N;
    }

    /** Append cell to the mesh */
    virtual void append(const SimpleCell<dim> &) {}

    /** Initialize markers with value */
    void init_markers(const int N, const int value) {
        require(N >= 0, "Invalid number of markers: " + d2s(N));
        markers = vector<int>(N, value);
    }
    
    /** Initialize marker appending */
    void init_markers(const int N) {
        require(N >= 0, "Invalid number of markers: " + d2s(N));
        markers.clear();
        markers.reserve(N);
    }

    /** Append cell marker */
    void append_marker(const int &m) {
        expect(get_n_markers() < static_cast<int>(markers.capacity()), "Allocated size of markers exceeded!");
        markers.push_back(m);
    }

    /** Set the number of available cells */
    void set_counter(const int n_cells) {
        require(n_cells >= 0, "Invalid number of cells: " + d2s(n_cells));
        i_cells = n_cells;
    }

    /** Copy the nodes from write buffer to read buffer */
    virtual void transfer(const bool write2read=false) {
        if(write2read)
            *n_cells_r = *n_cells_w;
        else
            *n_cells_w = *n_cells_r;
        i_cells = *n_cells_w;
    }

    /** Copy the cells from another mesh; mask can be used to copy only specified cells -
     * i-th true|false in mask means i-th cell is copied|skipped. */
    void copy(const TetgenCells<dim>& cells, const vector<bool>& mask={}) {
        const int n_cells = cells.size();

        // In case of empty or non-aligned mask, copy all the cells
        if (n_cells != static_cast<int>(mask.size())) {
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
        if (n_markers != static_cast<int>(mask.size())) {
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
        require(i >= 0 && i < static_cast<int>(get_n_markers()), "Invalid index: " + d2s(i));
        return markers[i];
    }

    /** Get cell centroid coordinates */
    Point3 get_centroid(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));

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
        require(node >= 0 && node < static_cast<int>(get_n_markers()), "Invalid index: " + d2s(node));
        markers[node] = m;
    }

    /** Calculate the neighbourlist for the nodes.
     * Two nodes are considered neighbours if they share a cell. */
    void calc_nborlist(vector<vector<unsigned>>& nborlist) const {
        nborlist = vector<vector<unsigned>>(get_n_nodes());
        for (SimpleCell<dim> cell : *this)
            for (int n1 : cell)
                for (int n2 : cell) {
                    if (n1 == n2) continue;
                    nborlist[n1].push_back(n2);
                }
    }

    /** Accessor for accessing i-th cell */
    SimpleCell<dim> operator [](const size_t i) const { return get_cell(i); }

    /** Iterator to access the cells */
    typedef Iterator<TetgenCells, SimpleCell<dim>> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    static constexpr int DIM = dim; //!< dimensionality of the cell; 1-node, 2-edge, 3-triangle etc

protected:
    int* n_cells_r;      ///< number of readable cells in mesh data
    int* n_cells_w;      ///< number of writable cells in mesh data
    tetgenio* reads;     ///< mesh data that has been processed by Tetgen
    tetgenio* writes;    ///< mesh data that will be fed to Tetgen

    int i_cells;         ///< cell counter
    vector<int> markers; ///< cell markers

    /** Return the cell type in vtk format;
     * 1-vertex, 3-line, 5-triangle, 9-quadrangle, 10-tetrahedron, 12-hexahedron */
    virtual int get_cell_type() const { return 0; }

    /** Return i-th cell */
    virtual SimpleCell<dim> get_cell(const int) const { return SimpleCell<dim>(); }

    /** Return number of readable nodes in the mesh */
    int get_n_nodes() const {
        return reads->numberofpoints;
    }

    /** Return i-th readable node from the mesh in point form */
    Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + d2s(i));
        const int n = n_coordinates * i;
        return Point3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
    }
    
    /** Return i-th readable node from the mesh in vector form */
    Vec3 get_vec(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + d2s(i));
        const int n = n_coordinates * i;
        return Vec3(reads->pointlist[n+0], reads->pointlist[n+1], reads->pointlist[n+2]);
    }

    /** Output mesh nodes and cells in .vtk format */
    void write_vtk_points_and_cells(ofstream &out) const {
        const size_t n_nodes = get_n_nodes();
        const size_t n_cells = size();
        const size_t celltype = get_cell_type();

        // Output the nodes
        out << "POINTS " << n_nodes << " double\n";
        for (size_t node = 0; node < n_nodes; ++node)
            out << get_node(node) << "\n";

        // Output the cells (tetrahedra, triangles, edges or vertices)
        out << "CELLS " << n_cells << " " << n_cells * (dim + 1) << "\n";
        for (size_t cl = 0; cl < n_cells; ++cl)
            out << dim << " " << get_cell(cl) << "\n";

        // Output cell types
        out << "CELL_TYPES " << n_cells << "\n";
        for (size_t cl = 0; cl < n_cells; ++cl)
            out << celltype << "\n";
    }

    /** Output mesh data in .vtk format */
    void write_vtk_data(ofstream &out) const {
        const size_t n_markers = get_n_markers();
        const size_t n_cells = size();

        // Output cell id-s
        out << "SCALARS ID int\nLOOKUP_TABLE default\n";
        for (size_t cl = 0; cl < n_cells; ++cl)
            out << cl << "\n";

        // Output cell markers
        if ((n_markers > 0) && (n_markers == n_cells)) {
            out << "SCALARS marker int\nLOOKUP_TABLE default\n";
            for (size_t cl = 0; cl < n_cells; ++cl)
                out << get_marker(cl) << "\n";
        }
    }

    /** Output mesh node data in .vtk format */
    void write_vtk_point_data(ofstream &out) const {
        if (get_cell_type() != VtkType::vertex) return;
        out << "POINT_DATA " << size() << "\n";
        write_vtk_data(out);
    }

    /** Output mesh cell data in .vtk format */
    void write_vtk_cell_data(ofstream &out) const {
        if (get_cell_type() == VtkType::vertex) return;
        out << "CELL_DATA " << size() << "\n";
        write_vtk_data(out);
    }

    /** Specify implemented output file formats */
    bool valid_extension(const string &ext) const {
        return ext == "vtk" || ext == "vtks";
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

    /** Modify the boundary value of i-th node */
    void set_boundary(const int i, const int value);

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
    void transfer(const bool write2read=true);

    /** Copy the nodes from another mesh */
    void copy(const TetgenNodes& nodes, const vector<bool>& mask={});

    /** Return the coordinates of i-th node as a 3D vector */
    Vec3 get_vec(const int i) const;

    /** Save the locations of the initially added nodes */
    void save_indices(const int n_surf, const int n_bulk, const int n_vacuum);

    /** Store the locations of different kinds of nodes that were produced while splitting tetrahedra into hexahedra */
    void save_hex_indices(const vector<int>& n_nodes);

    /** Calculate statistics about nodes */
    void calc_statistics();

    /** Transform nodes into Deal.II format */
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

        double xbox;    ///< distance between max & min x-coordinate
        double ybox;    ///< distance between max & min y-coordinate
        double zbox;    ///< distance between max & min z-coordinate
    } stat;

private:
    /** Return the vertex type in vtk format */
    int get_cell_type() const { return VtkType::vertex; }

    /** Return index of i-th node */
    SimpleCell<1> get_cell(const int i) const;

    void copy_statistics(const TetgenNodes& nodes);

    /** Write node data to .xyz file */
    void write_xyz(ofstream &out) const;

    /** Initialize statistics about nodes */
    void init_statistics();

    /** Specify implemented output file formats */
    bool valid_extension(const string &ext) const {
        return ext == "xyz" || ext == "movie" || ext == "vtk" || ext == "vtks";
    }
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

    /** Copy the nodes from one buffer to another */
    void transfer(const bool write2read=true);

    /** Delete the edges that are not on the perimeter of surface */
    void clean_sides(const Medium::Sizes& stat);

    /** Calculate statistics about ewrite_vtkdges */
    void calc_statistics();

    /** Struct holding statistics about edges */
    struct Stat {
        double edgemin;    //!< Minimum edge length
        double edgemax;    //!< Maximum edge length
    } stat;

private:
    /** Return the line type in vtk format */
    int get_cell_type() const { return VtkType::line; }

    /** Return i-th edge */
    SimpleCell<2> get_cell(const int i) const;

    /** Initialize statistics about edges */
    void init_statistics();
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

    /** Copy the nodes from one buffer to another */
    void transfer(const bool write2read=true);
    
    /** Copy the surface faces from another TetgenFaces */
    int copy_surface(const TetgenFaces& faces, const TetgenNodes::Stat& stat);

    /** Return the normal of i-th triangle */
    Vec3 get_norm(const int i) const;
    
    /** Return the area of i-th triangle */
    double get_area(const int i) const;

    /** Return the centroid of i-th triangle */
    Point3 get_centroid(const int i) const;

    /** Return indices of all tetrahedra that are connected to i-th triangle;
     * -1 means there's no tetrahedron */
    array<int,2> to_tets(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        const int I = 2 * i;
        return array<int,2>{reads->face2tetlist[I], reads->face2tetlist[I+1]};
    }

    /** Return indices of all quadrangles that are connected to i-th triangle*/
    array<int,3> to_quads(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        const int I = n_quads_per_tri * i;
        return array<int,3>{I, I+1, I+2};
    }

    /** Calculate statistics about triangles */
    void calc_statistics();

    /** Calculate the norms and areas for all the triangles */
    void calc_appendices();

    /** Struct holding statistics about triangles */
    struct Stat {
        double edgemin;    ///< min edge length
        double edgemax;    ///< max edge length
        double xmin;       ///< min x-coordinate of face centroids
        double xmax;       ///< max x-coordinate of face centroids
        double ymin;       ///< min y-coordinate of face centroids
        double ymax;       ///< max y-coordinate of face centroids
        double zmin;       ///< min z-coordinate of face centroids
        double zmax;       ///< max z-coordinate of face centroids
        double xbox;       ///< distance between max & min x-coordinate of face centroids
        double ybox;       ///< distance between max & min y-coordinate of face centroids
        double zbox;       ///< distance between max & min z-coordinate of face centroids
    } stat;

private:
    vector<double> areas;     ///< areas of triangles
    vector<Vec3> norms;       ///< norms of triangles
    vector<Point3> centroids; ///< pre-calculated centroids of triangles

    /** Return the triangle type in vtk format */
    int get_cell_type() const { return VtkType::triangle; }

    /** Return i-th face */
    SimpleCell<3> get_cell(const int i) const;

    /** Initialize statistics about triangles */
    void init_statistics();
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

    /** Copy the nodes from one buffer to another */
    void transfer(const bool write2read=true);

    /** Return indices of all triangles that are connected to i-th tetrahedron*/
    vector<int> to_tris(const int i) const {
        require(i >= 0 && i < (int)map2tris.size(), "Invalid index: " + d2s(i));
        return map2tris[i];
    }

    /** Return indices of all hexahedra that are connected to i-th tetrahedron*/
    array<int,4> to_hexs(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        const int I = n_hexs_per_tet * i;
        return array<int,4>{I, I+1, I+2, I+3};
    }

    /** Calculate statistics about tetrahedra */
    void calc_statistics();

    /** Store data for mapping tetrahedron to triangles */
    void store_map(const vector<vector<int>>& map) {
        map2tris = map;
    }

    /** Struct holding statistics about tetrahedra */
    struct Stat {
        double edgemin;    //!< Minimum edge length
        double edgemax;    //!< Maximum edge length
    } stat;

private:
    vector<vector<int>> map2tris; ///< data for mapping tetrahedron to the triangles

    /** Return the tetrahedron type in vtk format */
    int get_cell_type() const { return VtkType::tetrahedron; }

    /** Return i-th element */
    SimpleCell<4> get_cell(const int i) const;

    /** Initialize statistics about tetrahedra */
    void init_statistics();
};

/** Class for holding quadrangles that were generated from triangles */
class Quadrangles: public TetgenCells<4> {
public:
    Quadrangles() : TetgenCells<4>() {}
    Quadrangles(tetgenio *data) : TetgenCells<4> (data, &data->numberofvcells) {}
    Quadrangles(tetgenio *read, tetgenio *write) :
        TetgenCells<4>(read, write, &read->numberofvcells, &write->numberofvcells) {}

    /** Initialize quadrangles appending */
    void init(const int N);

    /** Append quadrangle to the mesh */
    void append(const SimpleCell<4> &cell);

    /** Get number of quadrangles in mesh */
    int size() const { return quads.size(); }

    /** Return indices of all hexahedra that are connected to i-th quadrangle;
     * -1 means there's no hexahedron */
    array<int,2> to_hexs(const int i) const {
        require(i >= 0 && i < (int)map2hexs.size(), "Invalid index: " + d2s(i));
        return map2hexs[i];
    }

    /** Return the index of triangle connected to i-th quadrangle */
    int to_tri(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return int(i / n_quads_per_tri);
    }

    /** Store data for mapping quadrangle to hexahedra */
    void store_map(const vector<array<int,2>>& map) {
        map2hexs = map;
    }

protected:
    vector<SimpleQuad> quads;
    vector<array<int,2>> map2hexs;

    /** Return the quadrangle type in vtk format */
    int get_cell_type() const { return VtkType::quadrangle; }

    /** Return i-th quadrangle */
    SimpleCell<4> get_cell(const int i) const;
};

/** Class for holding hexahedra that were generated from tetrahedra */
class Hexahedra: public TetgenCells<8> {
public:
    Hexahedra() : TetgenCells<8>() {}
    Hexahedra(tetgenio *data) : TetgenCells<8> (data, &data->numberofvcells) {}
    Hexahedra(tetgenio *read, tetgenio *write) :
        TetgenCells<8>(read, write, &read->numberofvcells, &write->numberofvcells) {}

    /** Initialize hexahedra appending */
    void init(const int N);

    /** Append hexahedron to the mesh */
    void append(const SimpleCell<8> &cell);

    /** Get number of hexahedra in mesh */
    int size() const { return hexs.size(); }

    /** Return the indices of quadrangles connected to i-th hexahedron */
    vector<int> to_quads(const int i) const {
        require(i >= 0 && i < (int)map2quads.size(), "Invalid index: " + d2s(i));
        return map2quads[i];
    }

    /** Return the index of tetrahedron connected to i-th hexahedron */
    int to_tet(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return int(i / n_hexs_per_tet);
    }

    /** Export vacuum hexahedra in Deal.II format */
    vector<dealii::CellData<3>> export_vacuum() const;

    /** Export bulk hexahedra in Deal.II format */
    vector<dealii::CellData<3>> export_bulk() const;

    /** Store the data for mapping hexahedron to quadrangles */
    void store_map(vector<vector<int>>& map) {
        map2quads = map;
    }

protected:
    vector<SimpleHex> hexs;
    vector<vector<int>> map2quads;

    /** Return the hexahedron type in vtk format */
    int get_cell_type() const { return VtkType::hexahedron; }

    /** Return i-th hexahedron */
    SimpleCell<8> get_cell(const int i) const;
};

} /* namespace femocs */

#endif /* TETGENCELLS_H_ */
