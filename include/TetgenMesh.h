/*
 * TetgenMesh.h
 *
 *  Created on: 3.10.2016
 *      Author: veske
 */

#ifndef TETGENMESH_H_
#define TETGENMESH_H_

#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "TetgenCells.h"
#include "Medium.h"

using namespace std;
namespace femocs {

/**
 * @brief Class to create and handle finite element mesh in Tetgen format,
 * http://wias-berlin.de/software/tetgen/1.5/
 *
 * The following routines use extensively Tetgen feature to append the new points to the
 * end of input point list while not modifying the ordering of input points.
 * This helps to distinguish surface nodes from the generated nodes in bulk and vacuum.
 */
class TetgenMesh {
public:
    TetgenMesh();
    ~TetgenMesh() {}

    /** Smoothen the triangles using different versions of Taubin smoothing algorithm */
    void smoothen(const int n_steps, const double lambda, const double mu, const string& algorithm);

    /** Function to generate simple mesh that consists of one element */
    int generate_simple();

    /** Function to generate mesh from surface, bulk and vacuum atoms */
    int generate(const Medium& bulk, const Medium& surf, const Medium& vacuum, const string& cmd);

    /** @brief Separate tetrahedra & triangles into hexahedra & quadrangles
     * Separation is done by adding node to the centroid of the tetrahedron edges,
     * faces and tetrahedron itself. These nodes are added to the end of input point list.
     * This helps to distinguish input "tetrahedral" and "triangular" nodes
     * from the appended "hexahedral" and "quadrangular" ones.
     */
    bool generate_hexahedra();
    
    /** Using the separated tetrahedra generate the triangular surface on the vacuum-material boundary */
    int generate_surface(const Medium::Sizes& sizes, const string& cmd1, const string& cmd2);

    /** Mark mesh nodes and elements by their location relative to the surface atoms */
    bool mark_mesh();

    /** Separate generated mesh into bulk and vacuum meshes */
    int separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd);

    /** Generate list of hexahedral nodes that surround the tetrahedral nodes. The resulting cells
     * resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells. */
    void calc_pseudo_3D_vorocells(vector<vector<unsigned>>& cells, const bool vacuum) const;

    /** Generate list of quadrangle nodes that surround the triangle nodes. The resulting cells
     * resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells. */
    void calc_pseudo_2D_vorocells(vector<vector<unsigned>>& cells) const;

    /** Use Tetgen built-in functions to write mesh elements data into vtk file */
    bool write(const string& file_name);

    /** Write bulk or vacuum mesh */
    void write_separate(const string& file_name, const int type);

    /** Delete the data of previously stored mesh and initialise a new one */
    void clear();

    /** Copy mesh from input to output or vice versa without modification */
    int transfer(const bool write2read=true);

    /** Perform Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd);

    /** Perform double Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd1, const string& cmd2);

    /** Map the triangle to the tetrahedron by specifying the region (vacuum or bulk)  */
    int tri2tet(const int tri, const int region) const;

    /** Map the quadrangle to the hexahedron by specifying the region (vacuum or bulk) */
    int quad2hex(const int quad, const int region) const;

    void test_mapping() const;

    TetgenNodes nodes = TetgenNodes(&tetIOout, &tetIOin); ///< data & operations for mesh nodes
    TetgenEdges edges = TetgenEdges(&tetIOout);           ///< data & operations for mesh edges
    TetgenFaces tris = TetgenFaces(&tetIOout, &tetIOin); ///< data & operations for mesh triangles
    TetgenElements tets = TetgenElements(&tetIOout, &tetIOin); ///< data & operations for mesh tetrahedra
    Quadrangles quads = Quadrangles(&tetIOout);           ///< data & operations for mesh quadrangles
    Hexahedra hexs = Hexahedra(&tetIOout);           ///< data & operations for mesh hexahedra

    static constexpr int n_coordinates = 3;     ///< # coordinates

    static constexpr int n_nodes_per_edge = 2; ///< # nodes on an edge
    static constexpr int n_nodes_per_tri = 3;  ///< # nodes on a triangle
    static constexpr int n_nodes_per_quad = 4; ///< # nodes on a quadrangle
    static constexpr int n_nodes_per_tet = 4;  ///< # nodes on a tetrahedron
    static constexpr int n_nodes_per_hex = 8;  ///< # nodes on a hexahedron

    static constexpr int n_edges_per_tri = 3;  ///< # edges on a triangle
    static constexpr int n_edges_per_quad = 4; ///< # edges on a quadrangle
    static constexpr int n_edges_per_tet = 6;  ///< # edges on a tetrahedron
    static constexpr int n_edges_per_hex = 12; ///< # edges on a hexahedron
    static constexpr int n_tris_per_tet = 4;   ///< # triangles on a tetrahedron
    static constexpr int n_tets_per_tri = 2;   ///< # tetrahedra connected to a triangle
    static constexpr int n_hexs_per_quad = 2;  ///< # hexahedra connected to a quadrangle
    static constexpr int n_hexs_per_tet = 4;   ///< # hexahedra connected to a tetrahedron
    static constexpr int n_quads_per_tri = 3;  ///< # quadrangles connected to a triangle
    static constexpr int n_quads_per_hex = 6;  ///< # quadrangles connected to a hexahedron

    /** String stream prints the statistics about the mesh */
    friend std::ostream& operator <<(std::ostream &s, const TetgenMesh &t) {
        s << "#hexs=" << t.hexs.size()
                << ", #tets=" << t.tets.size()
                << ", #quads=" << t.quads.size()
                << ", #tris=" << t.tris.size()
                << ", #edges=" << t.edges.size()
                << ", #nodes=" << t.nodes.size();
        return s;
    }

private:
    tetgenio tetIOin;   ///< Writable mesh data in Tetgen format
    tetgenio tetIOout;  ///< Readable mesh data in Tetgen format

    /** Smoothen the surface mesh using Taubin lambda|mu algorithm with inverse neighbour count weighting */
    void laplace_smooth(const double scale, const vector<vector<unsigned>>& nborlist);

    /** Smoothen the surface mesh using Taubin lambda|mu algorithm with Fujiwara weighting */
    void fujiwara_smooth(const double scale, const vector<vector<unsigned>>& nborlist);

    /** Smoothen the surface mesh using Taubin lambda|mu algorithm with curvature normal weighting */
    void curvature_norm_smooth(const double scale, const vector<vector<unsigned>>& nborlist);

    /** Locate the tetrahedron by the location of its nodes */
    int locate_element(SimpleElement& elem);

    /** @brief Group hexahedra & quadrangles around central tetrahedral & triangular node
     *
     * Each hexahedron has one and only one unique tetrahedral node, therefore
     * Hexahedra with the same tetrahedral node form the pseudo Voronoi cell of that node.
     * Similar thing is valid for quadrangles, just they are mapped to triangular nodes.
     * After the routine the sign of hexahedron marker will indicate whether the hex is in
     * vacuum (>0) or bulk (<0) and the (magnitude - 1) shows the index of its tetrahedral node.
     * The quadrangle marker shows the index of its triangular node.
     */
    void group_hexahedra();

    bool calc_ranks(vector<int>& ranks, const vector<vector<unsigned>>& nborlist);

    /** Calculate the mapping between quadrangle and hexahedron indices */
    void calc_quad2hex_mapping();

    /** Calculate the mapping between tetrahedron and triangle indices */
    void calc_tet2tri_mapping();

    /** Mark the tetrahedra by the location of nodes */
    void mark_elems();

    /** Mark the nodes by using DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
     * algorithm. The same algorithm is also used in cluster analysis. */
    bool mark_nodes();

    bool rank_and_mark_nodes();

    /** Mark the edges on the simulation cell perimeter by the node markers */
    void mark_edges();

    /** Mark the boundary faces of mesh */
    void mark_faces();

    /** Generate the edges from the elements.
     * Overlapping edges are not cleaned to make it easier to match them with element. */
    void generate_edges();

    /** Generate surface faces from elements and known location of surface nodes.
     * Overlapping faces are not cleaned for computational efficiency purposes. */
    void generate_manual_surface();
};

} /* namespace femocs */

#endif /* TETGENMESH_H_ */
