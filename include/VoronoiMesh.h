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
#include "Globals.h"
#include "FileWriter.h"

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

    /** Get number of nodes that are associated with the face */
    int size() const { return data->vfacetlist[id].elist[0]; }

    /** Get the area of the face.
     * See the theory in http://geomalgorithms.com/a01-_area.html#3D%20Polygons */
    double area();

    /** Get the centroid coordinates of the face */
    Point3 centroid();

    /** Return the neighbouring cell for the caller cell */
    int nborcell(const int caller_id) const;

    /** Accessor for accessing the index of i-th edge */
    int operator [](int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return data->vfacetlist[id].elist[i+1];
    }

    /** Iterator to access the indexes of edges */
    typedef Iterator<VoronoiFace, int> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    /** Stream for printing the nodes to file or console */
    friend std::ostream& operator <<(std::ostream &os, VoronoiFace &vf) {
        int node = -1;
        for (int edge : vf) {
            vf.get_node(edge, node);
            os << node << ' ';
        }

        return os;
    }

    /** Get the norm vector of the face */
    Vec3 norm(const int cell) const;

    /** Transform the node data from tetgenio into easily accessible form */
    void calc_verts();

    vector<Vec3> verts;  ///< coordinates of the face vertices

private:
    /** Calculate the unique node that is associated with the edge */
    void get_node(const int edge, int& node) const;
};

/** Class for accessing the Voronoi cell data */
class VoronoiCell : public Voronoi {
public:
    VoronoiCell() : Voronoi() {}
    VoronoiCell(tetgenio* data, const int i) : Voronoi(data, i) {}

    /** Get number of faces that make the cell */
    int size() const { return data->vcelllist[id][0]; }

    /** Accessor for accessing the i-th face */
    VoronoiFace operator [](int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return VoronoiFace(data, data->vcelllist[id][i+1]);
    }

    /** Get the indices of neighbouring Voronoi cells */
    vector<int> get_neighbours() const;

    /** Iterator to access the cell faces */
    typedef Iterator<VoronoiCell, VoronoiFace> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

private:
};

/** Virtual class for holding data that is common to Voronoi cells and Voronoi faces */
template<typename T>
class Voronois : public FileWriter {
public:
    Voronois() : tetio(NULL), _n_cells(NULL) {}
    Voronois(tetgenio* data, int* n_cells) : tetio(data), _n_cells(n_cells) {}
    virtual ~Voronois() {};

    /** Get number of cells in mesh */
    int size() const { return *_n_cells; }

    /** Initialize markers with default value */
    void init_markers() { markers = vector<int>(size(), TYPES.NONE); }

    /** Assign i-th marker */
    void set_marker(const int i, const int value) {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        markers[i] = value;
    }

    /** Return i-th marker */
    int get_marker(const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return markers[i];
    }

    /** Return pointer to markers */
    vector<int>* get_markers() { return &markers; }

    /** Accessor for accessing the i-th cell */
    T operator [](const int i) const {
        require(i >= 0 && i < size(), "Invalid index: " + d2s(i));
        return T(tetio, i);
    }

    /** Iterator to access the cells */
    typedef Iterator<Voronois, T> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

protected:
    tetgenio* tetio;     ///< mesh data that has been processed by Tetgen
    int* _n_cells;       ///< number of cells in mesh data
    vector<int> markers; ///< cell markers

    /** Return number of readable nodes in the mesh */
    int get_n_nodes() const { return tetio->numberofvpoints; }

    /** Return i-th node from the voronoi mesh */
    Point3 get_node(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + d2s(i));
        const int n = n_coordinates * i;
        return Point3(tetio->vpointlist[n+0], tetio->vpointlist[n+1], tetio->vpointlist[n+2]);
    }

    /** Return i-th node from the voronoi mesh in vector form */
    Vec3 get_vec(const int i) const {
        require(i >= 0 && i < get_n_nodes(), "Invalid index: " + d2s(i));
        const int n = n_coordinates * i;
        return Vec3(tetio->vpointlist[n+0], tetio->vpointlist[n+1], tetio->vpointlist[n+2]);
    }

    /** Specify file types that can be written */
    bool valid_extension(const string &ext) const {
        return ext == "vtk";
    }

    /** Write the header and point data of Voronois in .vtk format */
    void write_vtk(ofstream &out) const {
        // Output vtk file header
        FileWriter::write_vtk(out);

        // Output the nodes
        const size_t n_nodes = get_n_nodes();
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
    VoronoiFaces() : Voronois<VoronoiFace>() {}
    VoronoiFaces(tetgenio* data) : Voronois<VoronoiFace>(data, &data->numberofvfacets) {}

    /** Calculate the edge-neighbours - faces, that share the same edge */
    void calc_neighbours();

    /** Calculate and store centroids of all faces */
    void calc_centroids();

    /** Return the edge-neighbours of i-th Voronoi face */
    vector<int> get_neighbours(const int i) const {
        require(i >= 0 && i < static_cast<int>(neighbours.size()), "Invalid index: " + d2s(i));
        return neighbours[i];
    }

    /** Return the centroid of i-th Voronoi face */
    Vec3 get_centroid(const int i) const {
        require(i >= 0 && i < static_cast<int>(centroids.size()), "Invalid index: " + d2s(i));
        return centroids[i];
    }

private:
    vector<vector<int>> neighbours;  ///< edge-neighbours of faces
    vector<Vec3> centroids;          ///< centroids of faces

    /** Output in vtk format the Voronoi faces and the data associated with them */
    void write_cells(ofstream& out) const;
};

/** Class for accessing the Voronoi cells */
class VoronoiCells : public Voronois<VoronoiCell> {
public:
    VoronoiCells() : Voronois<VoronoiCell>() {}
    VoronoiCells(tetgenio* data) : Voronois<VoronoiCell>(data, &data->numberofvcells) {}

private:
    /** Output in vtk format the Voronoi cells and the data associated with them */
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
    int generate(const Medium& surface, const double latconst, const string& cmd1, const string& cmd2);

    /** Mark the cells and faces with nodes in the infinity */
    void clean();

    /** Extract the atoms and their areas whose Voronoi cells are exposed to vacuum */
    void extract_surface(Medium& surface, vector<Vec3>& areas, const Medium& nanotip);

    /** Calculate minimum z-coordinate of medium and mark mesh */
    void mark_mesh(const Medium& medium, const double latconst);

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

    /** Mark the Voronoi cells, faces and nodes by their relative location against top and bottom surfaces */
    void mark_mesh(const double zmin);

    /** Mark the cell and faces that are certainly on the surface */
    int mark_seed();

    /** Mark Voronoi faces that are on the vacuum-material boundary */
    void mark_faces(const double zmin, const int seed);

    /** Function to calculate face ranks to increase the accuracy of their marking */
    void calc_ranks(vector<int>& ranks, const int seedface);

    /** Mark Voronoi cells and nodes that are on the surface of material */
    void mark_cells_and_nodes();

    /** Perform triple Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd1, const string& cmd2, const string& cmd3);

    /** Perform double Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd1, const string& cmd2);

    /** Perform single Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd1);
};

} /* namespace femocs */

#endif /* VORONOIMESH_H_ */
