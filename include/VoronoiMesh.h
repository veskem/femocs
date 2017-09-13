/*
 * VoronoiMesh.h
 *
 *  Created on: 6.3.2017
 *      Author: veske
 */

#ifndef VORONOIMESH_H_
#define VORONOIMESH_H_

#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "Medium.h"
#include "TetgenCells.h"
#include "SolutionReader.h"
#include <fstream>

using namespace std;
namespace femocs {

/** Virtual class for holding data that is common to Voronoi cell and Voronoi face */
class Voronoi {
public:
    Voronoi() : id(-1), data(NULL) {}
    Voronoi(tetgenio* data, const int i) : id(i), data(data) {}
    virtual ~Voronoi() {}

    /** Get number of faces that make the cell */
    virtual int size() const { return 0; }

    int id;  ///< ID of Voronoi cell or face in mesh data list

protected:
    tetgenio* data;     ///< mesh data that has been processed by Tetgen
};

/** Class for accessing the Voronoi face data */
class VoronoiFace : public Voronoi {
public:
    VoronoiFace() : Voronoi() {}
    VoronoiFace(tetgenio* data, const int i) : Voronoi(data, i) {}
    ~VoronoiFace() {}

    /** Get number of nodes that are associated with the face */
    int size() const { return data->vfacetlist[id].elist[0]; }

    /** Get the area of the face in the form of its scaled norm.
     * See the theory in http://geomalgorithms.com/a01-_area.html#3D%20Polygons */
    Vec3 area();

    /** Get the centroid coordinates of the face */
    Vec3 centroid();

    /** Return the neighbouring cell for the caller cell */
    int nborcell(const int caller_id);

    /** Accessor for accessing the index of i-th edge */
    int operator [](int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return data->vfacetlist[id].elist[i+1];
    }

    /** Iterator to access the indexes of edges */
    typedef Iterator<VoronoiFace, int> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Stream for printing the nodes to file or console */
    friend std::ostream& operator <<(std::ostream &os, VoronoiFace &vf) {
        for (int edge : vf)
            os << vf.get_node(edge) << ' ';
        return os;
    }

    /** Calculate the unique node that is associated with the edge */
    int get_node(const int edge);

    /** Transform the node data from tetgenio into easily accessible form */
    void calc_verts();

    vector<Vec3> verts;  ///< coordinates of the face vertices

private:

    /** Get the norm vector of the face */
    Vec3 norm();
};

/** Class for accessing the Voronoi cell data */
class VoronoiCell : public Voronoi {
public:
    VoronoiCell() : Voronoi() {}
    VoronoiCell(tetgenio* data, const int i) : Voronoi(data, i) {}
    ~VoronoiCell() {}

    /** Get number of faces that make the cell */
    int size() const { return data->vcelllist[id][0]; }

    /** Accessor for accessing the i-th face */
    VoronoiFace operator [](int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return VoronoiFace(data, data->vcelllist[id][i+1]);
    }

    /** Get the indices of neighbouring Voronoi cells */
    vector<int> get_neighbours() const;

    Vec3 area() const;

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

    /** Get number of cells in mesh */
    int size() const { return *_n_cells; }

    /** Initialize markers with default value */
    void init_markers() { markers = vector<int>(size(), TYPES.NONE); }

    /** Assign i-th marker */
    void set_marker(const int i, const int value) {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        markers[i] = value;
    }

    /** Return i-th marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return markers[i];
    }

    /** Return pointer to markers */
    vector<int>* get_markers() { return &markers; }

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
    T operator [](const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return T(tetio, i);
    }

    /** Iterator to access the cells */
    typedef Iterator<Voronois, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

protected:
    static constexpr int n_coordinates = 3;  ///< number of spatial coordinates
    static constexpr int celltype = 7;       ///< vtk cell = polygon

    tetgenio* tetio;     ///< mesh data that has been processed by Tetgen
    int* _n_cells;       ///< number of cells in mesh data
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

    /** Write the header and point data of Voronois in .vtk format */
    void write_vtk(const string &file_name) const {
        const size_t n_nodes = get_n_nodes();

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

    /** Write the cell data (number and indices of vertices) to the file */
    virtual void write_cells(ofstream&) const {}
};

/** Class for accessing the faces of Voronoi cells */
class VoronoiFaces : public Voronois<VoronoiFace> {
public:
    /** VoronoiCells constructors */
    VoronoiFaces() : Voronois<VoronoiFace>() {}
    VoronoiFaces(tetgenio* data) : Voronois<VoronoiFace>(data, &data->numberofvfacets) {}

    void write_cells(ofstream& out) const;

    void calc_neighbours();

    void calc_centroids();

    vector<int> get_neighbours(const int i) const {
        require(i >= 0 && i < static_cast<int>(neighbours.size()), "Invalid index: " + to_string(i));
        return neighbours[i];
    }

    Vec3 get_centroid(const int i) const {
        require(i >= 0 && i < static_cast<int>(centroids.size()), "Invalid index: " + to_string(i));
        return centroids[i];
    }

private:
    vector<vector<int>> neighbours;
    vector<Vec3> centroids;
};

/** Class for accessing the Voronoi cells */
class VoronoiCells : public Voronois<VoronoiCell> {
public:
    /** VoronoiCells constructors */
    VoronoiCells() : Voronois<VoronoiCell>() {}
    VoronoiCells(tetgenio* data) : Voronois<VoronoiCell>(data, &data->numberofvcells) {}

    void write_cells(ofstream& out) const;
};

/** Class to calculate and handle Voronoi cells around (surface) atoms.
 * Voronoi cells are made with Tetgen, http://wias-berlin.de/software/tetgen/1.5/
 */
class VoronoiMesh {
public:
    VoronoiMesh();
    ~VoronoiMesh() {}

    /** Generate Voronoi cells around surface atoms */
    bool generate(const Medium& surface, const double latconst, const string& cmd1, const string& cmd2);

    bool generate_modi(const Medium& surface, const double latconst, const string& cmd1, const string& cmd2);

    /** Mark the cells and faces with nodes in the infinity */
    void clean();
    
    /** Extract the atoms and their areas whose Voronoi cells are exposed to vacuum */
    void extract_surface(Medium& surface, vector<Vec3>& areas, const Medium& nanotip);

    /** Calculate minimum z-coordinate of medium and mark mesh */
    void mark_mesh(const Medium& medium, const double latconst);

    /** Mark the Voronoi cells, faces and nodes by their relative location against top and bottom surfaces */
    void mark_mesh(const double zmin);

    /** Objects holding operations for accessing cell data */
    TetgenNodes nodes = TetgenNodes(&tetIOout, &tetIOin);
    TetgenElements elems = TetgenElements(&tetIOout, &tetIOin);
    VoronoiCells voros = VoronoiCells(&tetIOout);
    VoronoiFaces vfaces = VoronoiFaces(&tetIOout);

private:
    tetgenio tetIOin;   ///< Writable mesh data in Tetgen format
    tetgenio tetIOout;  ///< Readable mesh data in Tetgen format

    /** Find the Voronoi cell that for sure belongs to the surface */
    int get_seedcell();

    /** Mark the cell and faces that are certainly on the surface */
    int mark_seed();

    /** Mark Voronoi faces that are on the vacuum-material boundary */
    void mark_faces(const double zmin, const int seed);

    void mark_faces_vol2(const double zmin, const int seed);

    void calc_ranks(vector<int>& ranks, const int seedface);

    /** Mark Voronoi cells and nodes that are on the surface of material */
    void mark_cells_and_nodes();

    /** Perform triple Tetgen calculation on input buffer and store it in output one */
    bool recalc(const string& cmd1, const string& cmd2, const string& cmd3);

    bool recalc(const string& cmd1, const string& cmd2);

    bool recalc(const string& cmd1);
};

} /* namespace femocs */

#endif /* VORONOIMESH_H_ */
