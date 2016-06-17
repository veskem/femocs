/*
 * SolutionReader.h
 *
 *  Created on: 06.06.2016
 *      Author: veske
 */

#ifndef SOLUTIONREADER_H_
#define SOLUTIONREADER_H_

#include "Primitives.h"
#include "DealII.h"

using namespace std;
namespace femocs {

/** Class to extract solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader conctructor */
    SolutionReader();
    const void extract_solution(DealII* fem, Medium &medium);
    const void extract_statistics(Mesh &mesh);
    const void export_helmod(int n_atoms, double* Ex, double* Ey, double* Ez, double* Enorm);
    const void output(const string file_name);

private:
    DealII* fem;
    double longrange_efield;
    vector<Vec3> elfield;
    vector<double> elfield_norm;
    vector<double> potential;
    vector<double> face_qualities;
    vector<double> elem_qualities;

    /** Reserve memory for solution vectors */
    const void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);

    /** Map the indices of nodes to the indices of the elements and vertices
     * where the nodes are located in those elements.
     * The index of element and vertex are stored sequentally, so that
     * map[2*node_index] == element_index, map[2*node_index+1] == vertex_index
     */
    const vector<int> get_node2elem_map();

    const vector<int> get_node2vert_map();

    const vector<int> get_node2face_map(Mesh &mesh, int node);
    const vector<int> get_node2elem_map(Mesh &mesh, int node);

    /** Return mapping between Medium atoms and DealII mesh nodes.
     * Value -1 indicates that there's no node in the mesh that corresponds to the given atom.
     */
    const vector<int> get_medium2node_map(Medium &medium);
};

} /* namespace femocs */

#endif /* SOLUTIONREADER_H_ */
