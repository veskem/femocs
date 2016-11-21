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
#include "Medium.h"
#include "TetgenCells.h"

using namespace std;
namespace femocs {

/** Class to extract solution from DealII calculations */
class SolutionReader: public Medium {
public:
    /** SolutionReader conctructors */
    SolutionReader();

    /** Extract the electric potential and electric field values on the nodes of tetrahedra from FEM solution */
    const void extract_solution(DealII &fem, const TetgenNodes &nodes);

    const Solution get_solution(const int i) const;

    /** Electric field that is assigned to atoms not found from mesh.
     *  Its value is BIG to make it immediately visible from data set. */
    const double error_field = 1e20;

private:
    vector<Solution> solution;
    
    /** Reserve memory for solution vectors */
    const void reserve(const int n_nodes);

    /** Get i-th entry from all data vectors; i < 0 gives the header of data vectors */
    const string get_data_string(const int i);

    /** Return the mapping between tetrahedral and hexahedral meshes; -1 indicates that mapping for corresponding object was not found
     * @param fem         solution from Deal.II
     * @param tet2hex     mapping between tet- & hexmesh elements,
     * @param node2hex    mapping between tetmesh nodes & hexmesh elements,
     * @param node2vert   mapping between tetmesh nodes & hexmesh element's vertices. */
    const void get_maps(DealII& fem, vector<int>& tet2hex, vector<int>& node2hex, vector<int>& node2vert);
};

} // namespace femocs

#endif /* SOLUTIONREADER_H_ */
