/*
 * VoronoiCells.cpp
 *
 *  Created on: 18.1.2017
 *      Author: veske
 */

#include "VoronoiCells.h"

using namespace std;
namespace femocs  {

VoronoiCells::VoronoiCells() : nodes(NULL), hexahedra(NULL) {}

VoronoiCells::VoronoiCells(TetgenNodes* n, Hexahedra* h) : nodes(n), hexahedra(h) {}

// Initialize cell appending
void VoronoiCells::init(const int N) {
    cells.clear();
    cells.reserve(N);
}

// Append cell to the mesh
void VoronoiCells::append(const VoronoiCell& cell) {
    expect(get_n_cells() < cells.capacity(), "Allocated memory exceeded!");
    cells.push_back(cell);
}

// Initialize marker appending
void VoronoiCells::init_markers(const int N) {
    markers.clear();
    markers.reserve(N);
}

// Append cell marker
void VoronoiCells::append_marker(const int &m) {
    expect(get_n_markers() < markers.capacity(), "Allocated memory exceeded!");
    markers.push_back(m);
}

// Get number of cells in mesh
int VoronoiCells::get_n_cells() const {
    return cells.size();
}

// Get number of cell markers
int VoronoiCells::get_n_markers() const {
    return markers.size();
}

// Return i-th cell
VoronoiCell VoronoiCells::get_cell(const int i) const {
    require(i >= 0 && i < get_n_cells(), "Index out of bounds: " + to_string(i));
    return cells[i];
}

// Return i-th marker
int VoronoiCells::get_marker(const int i) const {
    require(i >= 0 && i < get_n_markers(), "Index out of bounds: " + to_string(i));
    return markers[i];
}

// Function to write cell data to file
void VoronoiCells::write(const string &file_name) const {

}


} /* namespace femocs */
