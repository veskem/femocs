/*
 * TetgenMesh.h
 *
 *  Created on: 3.10.2016
 *      Author: veske
 */

#ifndef TETGENMESH_H_
#define TETGENMESH_H_

#include "FileWriter.h"
#include "Macros.h"
#include "Primitives.h"
#include "Tetgen.h"
#include "TetgenCells.h"
#include "Medium.h"
#include "Config.h"

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
class TetgenMesh: public FileWriter {
public:
    TetgenMesh(const Config::Mesh* conf);
    ~TetgenMesh() {}

    /** Read mesh from file and generate mappings between cells */
    int read(const string &file, const string &cmd);

    /** Generate union mesh between generator points */
    int generate(const Medium& bulk, const Medium& surf, const Medium& vacuum);

    /** Generate simple mesh that consists of one element */
    int generate_simple();

    /** Separate generated mesh into bulk and vacuum meshes */
    int separate_meshes(TetgenMesh &bulk, TetgenMesh &vacuum, const string &cmd);

    /** Generate list of hexahedral nodes that surround the tetrahedral nodes. The resulting cells
     * resemble Voronoi cells but are still something else, i.e pseudo Voronoi cells. */
    void calc_pseudo_3D_vorocells(vector<vector<unsigned>>& cells, const bool vacuum) const;

    /** Map the triangle to the tetrahedron by specifying the region (vacuum or bulk)  */
    int tri2tet(const int tri, const int region) const;

    /** Map the quadrangle to the hexahedron by specifying the region (vacuum or bulk) */
    int quad2hex(const int quad, const int region) const;

    /** Write bulk or vacuum mesh */
    void write_separate(const string& file_name, const int type);

    TetgenNodes nodes = TetgenNodes(&tetIOout, &tetIOin); ///< data & operations for mesh nodes
    TetgenEdges edges = TetgenEdges(&tetIOout);           ///< data & operations for mesh edges
    TetgenFaces tris = TetgenFaces(&tetIOout, &tetIOin); ///< data & operations for mesh triangles
    TetgenElements tets = TetgenElements(&tetIOout, &tetIOin); ///< data & operations for mesh tetrahedra
    Quadrangles quads = Quadrangles(&tetIOout);           ///< data & operations for mesh quadrangles
    Hexahedra hexs = Hexahedra(&tetIOout);           ///< data & operations for mesh hexahedra

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

    /** Return mesh statistics as string */
    string to_str() const { stringstream ss; ss << (*this); return ss.str(); }

private:
    tetgenio tetIOin;   ///< Writable mesh data in Tetgen format
    tetgenio tetIOout;  ///< Readable mesh data in Tetgen format
    const Config::Mesh* conf;  ///< configuration data about mesh

    void write_bin(ofstream &out) const;

    void write_msh(ofstream &out) const;

    void write_vtk(ofstream &out);

    /** Function to generate mesh from surface, bulk and vacuum atoms */
    int generate_union(const Medium& bulk, const Medium& surf, const Medium& vacuum, const string& cmd);

    /** @brief Separate tetrahedra & triangles into hexahedra & quadrangles
     * Separation is done by adding node to the centroid of the tetrahedron edges,
     * faces and tetrahedron itself. These nodes are added to the end of input point list.
     * This helps to distinguish input "tetrahedral" and "triangular" nodes
     * from the appended "hexahedral" and "quadrangular" ones.
     */
    bool generate_hexahedra();

    /** Using the separated tetrahedra generate the triangular surface on the vacuum-material boundary */
    int generate_surface(const string& cmd1, const string& cmd2);

    /** Generate surface faces from elements and known location of surface nodes.
     * Overlapping faces are not cleaned. */
    void generate_manual_surface();

    /** Delete the data of previously stored mesh and initialise a new one */
    void clear();

    /** Copy mesh from input to output or vice versa without modification */
    int transfer(const bool write2read=true);

    /** Perform Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd);

    /** Perform double Tetgen calculation on input buffer and store it in output one */
    int recalc(const string& cmd1, const string& cmd2);

    /** Smoothen the triangles using different versions of Taubin smoothing algorithm */
    void smoothen();

    /** Smoothen the surface mesh using Taubin lambda|mu algorithm with inverse neighbour count weighting */
    void laplace_smooth(vector<Point3> &displacements, const double scale, const vector<vector<unsigned>>& nborlist);

    /** Smoothen the surface mesh using Taubin lambda|mu algorithm with Fujiwara weighting */
    void fujiwara_smooth(vector<Point3> &displacements, const double scale, const vector<vector<unsigned>>& nborlist);

    /** Handle the loss of mesh quality due to its smoothing */
    void restore_quality(vector<bool> &tet_map, vector<vector<unsigned>> &nborlist,
            const vector<Point3>& displacement1, const vector<Point3>&displacement2);

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

    /** @brief Calculate factors that show many connections the non-surface nodes have with its non-surface neighbors.
     *
     * Nodes with small amount of neighbours are either on the boundary of simubox or on the edge
     * of a hole, while nodes with large amount of neighbours are inside the bulk or vacuum domain.
     */
    bool calc_ranks(vector<int>& ranks, const vector<vector<unsigned>>& nborlist);

    /** Calculate the mapping between quadrangle and hexahedron indices */
    void calc_quad2hex2quad_mapping();

    /** Calculate the mapping between tetrahedron and triangle indices */
    int calc_tri2tet2tri_mapping(const string &cmd, int n_surf_faces);

    /** Read and calculate the mapping between tetrahedron and triangle indices */
    void read_tri2tet2tri_mapping(const string &filename);

    /** Calculate the mapping from tetrahedra to triangles */
    void calc_tet2tri_mapping();

    /** Mark the tetrahedra by the location of nodes */
    void mark_tets();

    /** Mark the boundary triangles of mesh */
    void mark_tris();

    /** Mark the edges on the simulation cell perimeter by the node markers */
    void mark_edges();

    /** @brief Mark the nodes by using the DBSCAN clustering algorithm
     *
     * DBSCAN (Density-Based Spatial Clustering of Applications with Noise) algorithm is
     * very sensitive against holes in the surface mesh that appear when there are
     * tetrahedra that have nodes both in vacuum and material domain.
     * To increase the tolerance against the holes in the surface, the ranking of the nodes
     * is calculated, i.e we measure how close the node is to the hole or to the simubox boundary.
     * Nodes with small ranking act as a boundary for the clustering algorithm.
     * The ranking helps to get rid of the effect of small holes.
     */
    bool rank_and_mark_nodes();

    /** Mark mesh nodes and elements by their location relative to the surface atoms */
    bool mark_mesh();

    /** Specify implemented output file formats */
    bool valid_extension(const string &ext) const {
        return ext == "restart" || ext == "bin" || ext == "msh";
    }
};

} /* namespace femocs */

#endif /* TETGENMESH_H_ */
