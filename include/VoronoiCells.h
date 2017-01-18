/*
 * VoronoiCells.h
 *
 *  Created on: 18.1.2017
 *      Author: veske
 */

#ifndef VORONOICELLS_H_
#define VORONOICELLS_H_


#include "Macros.h"
#include "Primitives.h"
#include "TetgenCells.h"

#include <fstream>

using namespace std;
namespace femocs {

class VoronoiCell {
public:
    VoronoiCell();

    void append(const unsigned int n);
    void append(const vector<unsigned int>& p);
    void init(const int n_nodes);

    /** Dimensionality of cell */
    int size() const { return nodes.size(); };

    /** Define the behaviour of string stream */
    friend std::ostream& operator <<(std::ostream &s, const VoronoiCell &t) {
        for (vector<unsigned int> polygon : t.polygons) {
            s << polygon.size() << ' ';
            for (unsigned int p : polygon) s << p << ' ';
            s << endl;
        }
        return s;
    }

    /** Return data as string */
    string to_str() const { stringstream ss; ss << this; return ss.str(); }

    /** Accessor for accessing the i-th node */
    const unsigned int& operator [](size_t i) const {
        require(i >= 0 && i < size(), "Invalid index: " + to_string(i));
        return nodes[i];
    }

    /** Iterator to access the cell nodes */
    typedef Iterator<VoronoiCell, unsigned int> iterator;
    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size()); }

    vector<vector<unsigned int>> polygons;
    vector<unsigned int> nodes;
};

class VoronoiCells {
public:
    /** Constructors */
    VoronoiCells();
    VoronoiCells(TetgenNodes* nodes, Hexahedra* hexahedra);

    /** Initialize cell appending */
    void init(const int N);

    /** Append cell to the mesh */
    void append(const VoronoiCell& cell);

    /** Initialize marker appending */
    void init_markers(const int N);

    /** Append cell marker */
    void append_marker(const int &m);

    /** Get number of cells in mesh */
    int get_n_cells() const;

    /** Get number of cell markers */
    int get_n_markers() const;

    /** Return i-th cell */
    VoronoiCell get_cell(const int i) const;

    /** Return i-th marker */
    int get_marker(const int i) const;

    /** Function to write cell data to file */
    void write(const string &file_name) const;

protected:
    TetgenNodes* nodes;
    Hexahedra* hexahedra;
    vector<int> markers;        ///< cell markers
    vector<VoronoiCell> cells;  ///< voronoi cells

    /** Output cells in .vtk format */
    void write_vtk(const string &file_name) const;
};

} /* namespace femocs */

#endif /* VORONOICELLS_H_ */
