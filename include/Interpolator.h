/*
 * Interpolator.h
 *
 *  Created on: 19.10.2016
 *      Author: veske
 */

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_


#include "Primitives.h"
#include "TetgenMesh.h"
#include "InterpolatorCells.h"
#include "CurrentsAndHeating.h"
#include "CurrentsAndHeatingStationary.h"
#include "Laplace.h"

using namespace std;
namespace femocs {

/** General class for interpolating solution inside mesh
 * Class holds subclasses (InterpolatorCells) that hold the actual interpolation data and routines
 * that allow performing interpolation either on surface or in space.
 * 
 * For tetrahedral and triangular interpolation, barycentric coordinates (BCC-s)
 * are used to calculate the shape (aka interpolation) functions.
 * Compact theory how to find BCC-s:
 * http://steve.hollasch.net/cgindex/geometry/ptintet.html
 *
 * Properties of determinant:
 *   http://www.vitutor.com/alg/determinants/properties_determinants.html
 *   http://www.vitutor.com/alg/determinants/minor_cofactor.html
 *
 * c++ code to find and handle BCC-s:
 *   http://dennis2society.de/painless-tetrahedral-barycentric-mapping
 *
 * Theory about shape functions inside different primitives.
 * Linear & quadratic triangle & quadrangle:
 *   https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch18.d/IFEM.Ch18.pdf
 * Linear tetrahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
 * Quadratic tetrahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
 * Linear hexahedron:
 *   https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
 */
class Interpolator {
public:
    Interpolator(const string& norm_label, const string& scalar_label);
    ~Interpolator() {};

    /** Initialise interpolator and store solution with default value */
    void initialize(const TetgenMesh* mesh, const double empty_value=0);

    /** Extract the electric potential and field values from FEM solution */
    void extract_solution(fch::Laplace<3>* fem);

    /** Extract the current density and stationary temperature values from FEM solution */
    void extract_solution(fch::CurrentsAndHeatingStationary<3>* fem);

    /** Extract the current density and transient temperature values from FEM solution */
    void extract_solution(fch::CurrentsAndHeating<3>& fem);

    int update_point_cell(Point3& point, int current_cell, bool deal_index = true);

    InterpolatorNodes nodes;      ///< vertices and solutions on them
    LinearTriangles lintris;      ///< data & operations for linear triangular interpolation
    LinearTetrahedra lintets;     ///< data & operations for linear tetrahedral interpolation
    QuadraticTriangles quadtris;  ///< data & operations for quadratic triangular interpolation
    QuadraticTetrahedra quadtets; ///< data & operations for quadratic tetrahedral interpolation
    LinearQuadrangles linquads;   ///< data & operations for linear quadrangular interpolation
    LinearHexahedra linhexs;      ///< data & operations for linear hexahedral interpolation

private:
    const TetgenMesh* mesh;         ///< Full mesh data with nodes, faces, elements etc
    int empty_value;                ///< Solution value for nodes outside the Deal.II mesh

    /** Calculate the mapping between Femocs & deal.II mesh nodes,
     *  nodes & hexahedral elements and nodes & element's vertices.
     *  -1 indicates that mapping for corresponding object was not found */
    void get_maps(vector<int>& femocs2deal, vector<int>& cell_indxs, vector<int>& vert_indxs,
            dealii::Triangulation<3>* tria, dealii::DoFHandler<3>* dofh) const;

    /** Transfer solution from FEM solver to Interpolator */
    void store_solution(const vector<int>& femocs2deal,
            const vector<dealii::Tensor<1, 3>> vec_data, const vector<double> scal_data);

    /** Force the solution on tetrahedral nodes to be the weighed average of the solutions on its
     *  surrounding hexahedral nodes */
    bool average_sharp_nodes(const bool vacuum);
};

} // namespace femocs

#endif /* INTERPOLATOR_H_ */
